"""
Piuma / MTS Contact-Model Data Analyzer - Full Loading Hertz Fit Version
Date: 2026-08-25

Main updates in this version:
1. Hertz fitting uses the entire first loading branch, from recording start to first prominent peak.
2. Sliding-window selection is removed to avoid bias toward a small local region.
3. The contact-boundary x0 and force-offset C bounds are wide, and multi-start fitting is used.
4. The model supports either crosshead direction:
       F = C + (4/3)E*sqrt(R)*max(0, direction*(x - x0))^(3/2)
5. The output table includes Hertz_fit_R2, Hertz_fit_RMSE_uN, and Hertz_fit_NRMSE_percent.
6. The plot shades the full loading branch used for fitting and overlays the Hertz fit over that branch.
7. The GUI has an Exclude/Include Selected button; excluded rows are marked in the table and omitted from Excel export.

Units:
- Crosshead: nm, converted to m for fitting
- Probe/sample radius: nm, converted to m
- Time: s
- Load: uN, converted to N for fitting
- E output: Pa
"""

from tkinter import filedialog
import os
import numpy as np
import pickle
import tkinter as tk
from tkinter import ttk, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from scipy.signal import find_peaks


# ---------------------------- Models & Utilities ----------------------------

def hertzian_contact_signed(x_m, E_Pa, R_m, x0_m, direction=1.0, C_N=0.0):
    """
    Hertz contact model with direction and vertical force offset.

    direction = +1:
        indentation increases as x increases.

    direction = -1:
        indentation increases as x decreases.

    F = C + (4/3) E sqrt(R) max(0, direction*(x - x0))^(3/2)
    """
    indentation_m = np.maximum(0.0, direction * (x_m - x0_m))
    return C_N + (4.0 / 3.0) * E_Pa * np.sqrt(R_m) * indentation_m ** 1.5


def bilateral_filter_1d(signal, radius=10, sigma_r=0.05, sigma_d=5.0):
    signal = np.asarray(signal, dtype=float)
    out = np.zeros_like(signal)
    n = len(signal)

    for i in range(n):
        weight_sum = 0.0
        acc = 0.0

        for k in range(-radius, radius + 1):
            j = min(max(i + k, 0), n - 1)
            range_weight = np.exp(-((signal[j] - signal[i]) ** 2) / (2 * sigma_r ** 2))
            spatial_weight = np.exp(-(k * k) / (2 * sigma_d ** 2))
            weight = range_weight * spatial_weight
            weight_sum += weight
            acc += signal[j] * weight

        out[i] = acc / (weight_sum if weight_sum else 1.0)

    return out


def read_file(folder_path):
    """
    Read each .txt file in the first subfolder, if present.

    Returns list of dictionaries with:
        filename, times, load, crosshead
    """
    recs = []

    if not os.path.isdir(folder_path):
        return recs

    subs = [e for e in os.listdir(folder_path) if os.path.isdir(os.path.join(folder_path, e)) and e.lower() != "calibration"]

    if subs:
        folder_path = os.path.join(folder_path, subs[0])

    for filename in os.listdir(folder_path):
        if not filename.lower().endswith(".txt") or os.path.splitext(filename)[0].lower() == "calibration":
            continue

        fp = os.path.join(folder_path, filename)
        c, l, t = [], [], []
        read_flag = False

        try:
            with open(fp, "r", encoding="utf-8-sig", errors="ignore") as f:
                for raw in f:
                    line = raw.strip()

                    if not read_flag and line.startswith("Time") and "\t" in line:
                        read_flag = True
                        continue

                    if read_flag and line:
                        parts = line.split("\t")

                        if len(parts) >= 3:
                            try:
                                t.append(float(parts[0]))
                                l.append(float(parts[1]))
                                c.append(float(parts[2]))
                            except ValueError:
                                pass

            if t and l and c:
                recs.append({
                    "filename": filename[:-4],
                    "times": t,
                    "load": l,
                    "crosshead": c,
                })

        except Exception as e:
            print("Read failed:", fp, e)

    return recs


# ------------------------------ Data Container ------------------------------

class Data:
    """
    One indentation repetition.

    Inputs:
        crosshead_nm: nm
        load_uN: uN
        times_s: s
        radius_nm: nm

    Internally:
        crosshead -> m
        load -> N
        radius -> m
    """

    def __init__(self, crosshead_nm, load_uN, times_s, radius_nm, name, out_dir):
        self.x_nm = np.asarray(crosshead_nm, dtype=float)

        # Preserve your previous behavior: shift the initial load upward if it starts negative.
        # This avoids negative loads in the Hertz fit, but still plots the same shifted force.
        load_uN = np.asarray(load_uN, dtype=float)
        self.f_uN = load_uN + abs(load_uN[0])

        self.t_s = np.asarray(times_s, dtype=float)
        self.R_nm = float(radius_nm)
        self.name = name
        self.out_dir = out_dir

        self.x_m = self.x_nm * 1e-9
        self.f_N = self.f_uN * 1e-6
        self.R_m = self.R_nm * 1e-9

        self.filt_uN = None

    # ---------- Helpers ----------

    def first_peak_index(self, min_prom_frac=0.05, min_distance_frac=0.01):
        """Return the index of the first prominent load peak."""
        y = np.asarray(self.f_uN, dtype=float)
        n = len(y)

        if n == 0:
            return 0

        y_range = np.nanmax(y) - np.nanmin(y)
        prominence = max(1e-12, min_prom_frac * y_range)
        distance = max(1, int(n * min_distance_frac))

        peaks, _ = find_peaks(y, prominence=prominence, distance=distance)

        return int(peaks[0]) if peaks.size else int(np.nanargmax(y))

    def filter_load(self):
        if self.filt_uN is None:
            self.filt_uN = bilateral_filter_1d(self.f_uN, radius=10, sigma_r=0.05, sigma_d=5.0)

    def plasticity_points(self):
        """
        Baseline in [10,50] s, threshold = baseline + 0.008 uN.
        Find crossings from below to above in [60,65] s.

        Returns:
            idx1, idx2, plasticity_nm
        """
        t = self.t_s
        f = self.f_uN

        if len(t) == 0:
            return None, None, np.nan

        mask = (t >= 10) & (t <= 50)
        baseline = np.mean(f[mask]) if np.any(mask) else 0.0
        threshold = baseline + 0.008

        win = (t >= 60) & (t <= 65)
        idx = np.where((f[:-1] < threshold) & (f[1:] >= threshold) & win[1:])[0]

        if idx.size == 0:
            return 0, len(t) - 1, np.nan

        zero1 = 0
        zero2 = int(idx.max())
        plast_nm = self.x_nm[zero2] - self.x_nm[zero1]

        return zero1, zero2, plast_nm

    @staticmethod
    def _fit_metrics(y_true, y_pred):
        y_true = np.asarray(y_true, dtype=float)
        y_pred = np.asarray(y_pred, dtype=float)

        residual = y_true - y_pred
        ss_res = float(np.sum(residual ** 2))
        ss_tot = float(np.sum((y_true - np.mean(y_true)) ** 2))

        r2 = np.nan if ss_tot <= 0 else 1.0 - ss_res / ss_tot
        rmse = float(np.sqrt(np.mean(residual ** 2)))

        y_range = float(np.nanmax(y_true) - np.nanmin(y_true))
        nrmse_pct = np.nan if y_range <= 0 else 100.0 * rmse / y_range

        return r2, rmse, nrmse_pct

    @staticmethod
    def _safe_corr(x, y):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        ok = np.isfinite(x) & np.isfinite(y)

        if np.sum(ok) < 3:
            return np.nan

        x = x[ok]
        y = y[ok]

        if np.nanstd(x) == 0 or np.nanstd(y) == 0:
            return np.nan

        return float(np.corrcoef(x, y)[0, 1])

    # ---------- Full loading-branch Hertz fit ----------

    def hertz_fit_linearized_loading_with_C(
        self,
        w_min=45,
        w_max_frac=1.0,
        r2_floor=0.0,
        force_low_frac=0.0,
        force_high_frac=1.0,
        x0_bound_span_factor=10.0,
        C_bound_span_factor=10.0,
    ):
        """
        Fit the ENTIRE first loading branch with the Hertz contact model.

        This version intentionally does NOT use a sliding/local window, because
        local-window selection can bias the result toward a small region of the
        loading curve.

        Fit region:
            from the beginning of the recording to the first prominent load peak.

        Model:
            F = C + (4/3) E sqrt(R) max(0, direction * (x - x0))^(3/2)

        Boundary changes:
            x0 and C are given wide bounds so the optimizer has enough room to
            place the contact point and vertical force offset.

        Returns:
            E_Pa,
            popt = (E_Pa, x0_m, C_N, direction),
            fit_idx = (i0, i1),
            fit_r2,
            fit_rmse_uN,
            fit_nrmse_pct
        """
        n = len(self.x_m)

        if n < 12 or self.R_m <= 0:
            return np.nan, None, (0, n), np.nan, np.nan, np.nan

        # 1. Use the full first loading branch: start -> first peak.
        ipk = self.first_peak_index()
        ipk = int(np.clip(ipk, 8, n - 1))

        i0_orig = 0
        i1_orig = ipk + 1

        x_fit = self.x_m[i0_orig:i1_orig]
        f_fit = self.f_N[i0_orig:i1_orig]

        valid = np.isfinite(x_fit) & np.isfinite(f_fit)
        x_seg = x_fit[valid]
        f_seg = f_fit[valid]

        if len(x_seg) < 8:
            return np.nan, None, (i0_orig, i1_orig), np.nan, np.nan, np.nan

        # 2. Determine crosshead direction over the full loading branch.
        corr = self._safe_corr(x_seg, f_seg)

        if np.isfinite(corr):
            direction = 1.0 if corr >= 0 else -1.0
        else:
            dx = x_seg[-1] - x_seg[0]
            df = f_seg[-1] - f_seg[0]
            direction = 1.0 if np.sign(dx) == np.sign(df) else -1.0

        x_span = float(np.nanmax(x_seg) - np.nanmin(x_seg))
        f_span = float(np.nanmax(f_seg) - np.nanmin(f_seg))

        if x_span <= 0 or f_span <= 0:
            return np.nan, None, (i0_orig, i1_orig), np.nan, np.nan, np.nan

        indent_span = max(1e-12, abs(direction * (x_seg[-1] - x_seg[0])))

        # 3. Broad initial modulus guess.
        E_guess = f_span / (
            (4.0 / 3.0) * np.sqrt(self.R_m) * indent_span ** 1.5 + 1e-30
        )

        if not np.isfinite(E_guess) or E_guess <= 0:
            E_guess = 1_000.0

        # 4. Wide bounds for x0 and C.
        x_min = float(np.nanmin(x_seg))
        x_max = float(np.nanmax(x_seg))
        f_min = float(np.nanmin(f_seg))
        f_max = float(np.nanmax(f_seg))

        x0_lb = x_min - x0_bound_span_factor * x_span
        x0_ub = x_max + x0_bound_span_factor * x_span

        C_lb = f_min - C_bound_span_factor * f_span
        C_ub = f_max + C_bound_span_factor * f_span

        lower_bounds = np.array([1e-12, x0_lb, C_lb], dtype=float)
        upper_bounds = np.array([np.inf, x0_ub, C_ub], dtype=float)

        # 5. Multi-start fitting. Still fits the full loading curve every time.
        #    This helps avoid local minima without reintroducing sliding-window bias.
        if direction > 0:
            x0_base = x_seg[0]
        else:
            x0_base = x_seg[0]

        x0_guesses = [
            x0_base,
            x_min,
            x_max,
            x_min - 0.25 * x_span,
            x_max + 0.25 * x_span,
            x_seg[0] - direction * 0.25 * x_span,
            x_seg[0] - direction * 1.00 * x_span,
        ]

        C_guesses = [
            f_seg[0],
            f_min,
            0.0,
            np.median(f_seg[:max(3, len(f_seg)//10)]),
        ]

        E_guesses = [
            E_guess * 0.01,
            E_guess * 0.1,
            E_guess,
            E_guess * 10.0,
            E_guess * 100.0,
            1_000.0,
            10_000.0,
        ]

        def residuals(p):
            E, x0, C = p
            return hertzian_contact_signed(
                x_seg,
                E,
                self.R_m,
                x0,
                direction=direction,
                C_N=C,
            ) - f_seg

        best = None

        for Eg in E_guesses:
            for x0g in x0_guesses:
                for Cg in C_guesses:
                    p0 = np.array([Eg, x0g, Cg], dtype=float)

                    if not np.all(np.isfinite(p0)):
                        continue

                    p0[0] = max(p0[0], lower_bounds[0] * 10)
                    p0[1] = np.clip(p0[1], lower_bounds[1] + 1e-18, upper_bounds[1] - 1e-18)
                    p0[2] = np.clip(p0[2], lower_bounds[2] + 1e-18, upper_bounds[2] - 1e-18)

                    try:
                        res = least_squares(
                            residuals,
                            x0=p0,
                            bounds=(lower_bounds, upper_bounds),
                            loss="linear",
                            max_nfev=80000,
                        )
                    except Exception:
                        continue

                    if not res.success:
                        continue

                    E_fit, x0_fit, C_fit = res.x

                    f_pred = hertzian_contact_signed(
                        x_seg,
                        E_fit,
                        self.R_m,
                        x0_fit,
                        direction=direction,
                        C_N=C_fit,
                    )

                    fit_r2, rmse_N, nrmse_pct = self._fit_metrics(f_seg, f_pred)

                    if not np.isfinite(rmse_N):
                        continue

                    score = rmse_N

                    if best is None or score < best["score"]:
                        best = {
                            "score": float(score),
                            "E": float(E_fit),
                            "x0": float(x0_fit),
                            "C": float(C_fit),
                            "direction": float(direction),
                            "r2": float(fit_r2) if np.isfinite(fit_r2) else np.nan,
                            "rmse_N": float(rmse_N),
                            "nrmse_pct": float(nrmse_pct) if np.isfinite(nrmse_pct) else np.nan,
                        }

        if best is None:
            return np.nan, None, (i0_orig, i1_orig), np.nan, np.nan, np.nan

        popt = (
            best["E"],
            best["x0"],
            best["C"],
            best["direction"],
        )

        return (
            best["E"],
            popt,
            (i0_orig, i1_orig),
            best["r2"],
            best["rmse_N"] * 1e6,
            best["nrmse_pct"],
        )

    # ---------- Plotting ----------

    def _save_fig_pickle(self, fig, subdir):
        folder = os.path.join(self.out_dir, "Figure", subdir)
        os.makedirs(folder, exist_ok=True)
        path = os.path.join(folder, f"{self.name}.pkl")

        with open(path, "wb") as fh:
            pickle.dump(fig, fh)

        plt.close(fig)

    def plot_time_vs_load(self):
        self.filter_load()
        fig, ax = plt.subplots()
        ttl = f"{self.name} time vs load"
        ax.plot(self.t_s, self.filt_uN, label="Filtered Load")
        ax.set_title(ttl)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Load (uN)")
        ax.legend(loc="best")
        self._save_fig_pickle(fig, "time VS load")

    def plot_time_vs_crosshead(self):
        fig, ax = plt.subplots()
        ttl = f"{self.name} time vs crosshead"
        ax.plot(self.t_s, self.x_nm, label="Crosshead")
        ax.set_title(ttl)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Crosshead (nm)")
        ax.legend(loc="best")
        self._save_fig_pickle(fig, "time VS crosshead")

    def plot_crosshead_vs_load_full_with_segment_fit(
        self,
        popt=None,
        fit_idx=(0, None),
        plast_pts=None,
        fit_r2=np.nan,
        fit_rmse_uN=np.nan,
        fit_nrmse_pct=np.nan,
    ):
        """
        Plot full raw curve and overlay Hertz fit over the full loading branch used for fitting.
        """
        fig, ax = plt.subplots()
        ttl = f"{self.name} crosshead vs load"

        ax.plot(self.x_nm, self.f_uN, label="Raw")

        # Plasticity markers.
        if plast_pts and all(p is not None for p in plast_pts[:2]):
            p1, p2 = int(plast_pts[0]), int(plast_pts[1])

            if 0 <= p1 < len(self.x_nm):
                ax.plot(self.x_nm[p1], self.f_uN[p1], "bo", label="Plasticity start")

            if 0 <= p2 < len(self.x_nm):
                ax.plot(self.x_nm[p2], self.f_uN[p2], "ro", label="Plasticity second")

        # Segment-only Hertz fit overlay.
        if popt is not None and fit_idx is not None:
            i0, i1 = fit_idx
            i1 = len(self.x_nm) if i1 is None else i1
            i0 = max(0, int(i0))
            i1 = min(len(self.x_nm), int(i1))

            if i1 > i0 + 2:
                E, x0_m, C_N, direction = popt
                x_seg_nm = self.x_nm[i0:i1]

                y_fit_N = hertzian_contact_signed(
                    self.x_m[i0:i1],
                    E,
                    self.R_m,
                    x0_m,
                    direction=direction,
                    C_N=C_N,
                )

                ax.axvspan(
                    min(self.x_nm[i0], self.x_nm[i1 - 1]),
                    max(self.x_nm[i0], self.x_nm[i1 - 1]),
                    alpha=0.08,
                    label="Fit window",
                )

                label = "Hertz fit"
                if np.isfinite(fit_r2):
                    label += f" (R²={fit_r2:.3f})"

                ax.plot(x_seg_nm, y_fit_N * 1e6, linewidth=2.5, label=label)

                info_text = f"E = {E:.4g} Pa"
                info_text += f"\nR² = {fit_r2:.4g}" if np.isfinite(fit_r2) else "\nR² = NA"

                if np.isfinite(fit_rmse_uN):
                    info_text += f"\nRMSE = {fit_rmse_uN:.3g} µN"

                if np.isfinite(fit_nrmse_pct):
                    info_text += f"\nNRMSE = {fit_nrmse_pct:.2f}%"

                info_text += f"\nDir = {direction:+.0f}"

                ax.text(
                    0.02,
                    0.98,
                    info_text,
                    transform=ax.transAxes,
                    verticalalignment="top",
                    bbox=dict(boxstyle="round", alpha=0.2),
                )

        ax.set_title(ttl)
        ax.set_xlabel("Crosshead (nm)")
        ax.set_ylabel("Load (uN)")
        ax.legend(loc="best")
        self._save_fig_pickle(fig, "crosshead VS load")


# ------------------------------ Controller/GUI ------------------------------

class Controller:
    def __init__(self, data_root):
        self.data_root = data_root
        self.pkl_dir = None
        self.test_names = []
        self.current_index = 0
        self.tree = None

    def plot_current(self):
        if not self.pkl_dir or not self.test_names:
            return

        name = self.test_names[self.current_index]
        pkl = os.path.join(self.data_root, "Figure", self.pkl_dir, f"{name}.pkl")

        for w in left_frame.winfo_children():
            w.destroy()

        if not os.path.isfile(pkl):
            tk.Label(left_frame, text=f"No figure for {name}").pack()
            return

        try:
            with open(pkl, "rb") as fh:
                fig = pickle.load(fh)
        except Exception as e:
            tk.Label(left_frame, text=f"Load failed for {name}: {e}").pack()
            return

        canvas = FigureCanvasTkAgg(fig, master=left_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        NavigationToolbar2Tk(canvas, left_frame).update()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    def set_plot_dir(self, subfolder):
        self.pkl_dir = subfolder
        self.plot_current()

    def next_test(self):
        if not self.test_names:
            return
        self.current_index = (self.current_index + 1) % len(self.test_names)
        self.plot_current()

    def prev_test(self):
        if not self.test_names:
            return
        self.current_index = (self.current_index - 1) % len(self.test_names)
        self.plot_current()

    def on_row_double_click(self, _event):
        sel = self.tree.selection()
        if not sel:
            return

        vals = self.tree.item(sel[0], "values")
        if not vals:
            return

        name = vals[0]
        if name in self.test_names:
            self.current_index = self.test_names.index(name)
            self.plot_current()


# ------------------------------- App Functions ------------------------------

def create_popup(text="Export Successful"):
    popup = tk.Toplevel(root)
    popup.title("Notification")
    popup.geometry("420x140")
    tk.Label(popup, text=text, font=("Helvetica", 12), wraplength=380).pack(pady=10)
    tk.Button(popup, text="Close", command=popup.destroy).pack(pady=5)


# Manual overrides are kept separate so applying a new cutoff never undoes
# the user's explicit Exclude/Include Selected choices.
manual_exclusion_overrides = {}


def get_r2_cutoff():
    """Get optional R² minimum for Excel export (blank means no cutoff)."""
    raw = entry_r2_cutoff.get().strip()
    if not raw:
        return None
    try:
        cutoff = float(raw)
    except ValueError as exc:
        raise ValueError("Enter a numeric R² cutoff between 0 and 1, or leave it blank.") from exc
    if not np.isfinite(cutoff) or not (0.0 <= cutoff <= 1.0):
        raise ValueError("R² cutoff must be between 0 and 1, or left blank.")
    return cutoff


def apply_r2_cutoff(refresh_table=True):
    """Exclude fits below the cutoff (or missing R²); preserve manual choices."""
    global result_df
    try:
        cutoff = get_r2_cutoff()
    except ValueError as exc:
        messagebox.showwarning("R² cutoff", str(exc))
        return False

    if result_df is None or result_df.empty:
        return True

    auto_excluded = pd.Series(False, index=result_df.index, dtype=bool)
    if cutoff is not None:
        r2 = pd.to_numeric(result_df["Hertz_fit_R2"], errors="coerce")
        auto_excluded = (r2 < cutoff) | ~np.isfinite(r2)

    for idx, manual_choice in manual_exclusion_overrides.items():
        if idx in auto_excluded.index:
            auto_excluded.loc[idx] = manual_choice

    result_df["Exclude_from_Excel"] = auto_excluded.astype(bool)
    if refresh_table and controller is not None:
        update_table(result_df, controller)
    return True


def export_to_excel():
    global result_df

    if result_df is None or result_df.empty:
        messagebox.showwarning("Export", "No data to export.")
        return

    # Use the current cutoff even if it was edited after pressing Run,
    # and refresh the table so the highlighted exclusions match the export.
    if not apply_r2_cutoff(refresh_table=True):
        return

    out_path = os.path.join(selected_dir_var.get() or ".", "output.xlsx")

    try:
        export_df = result_df.copy()

        if "Exclude_from_Excel" in export_df.columns:
            keep_mask = ~export_df["Exclude_from_Excel"].astype(bool)
            excluded_n = int((~keep_mask).sum())
            export_df = export_df.loc[keep_mask].copy()
            export_df = export_df.drop(columns=["Exclude_from_Excel"])
        else:
            excluded_n = 0

        with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
            export_df.to_excel(writer, sheet_name="Results", index=False)

        create_popup(
            f"Exported to:\n{out_path}\n\n"
            f"Excluded rows not exported: {excluded_n}"
        )
    except Exception as e:
        messagebox.showerror("Export failed", str(e))


def toggle_exclude_selected():
    """
    Toggle whether the selected table row should be excluded from Excel export.

    The row stays visible in the GUI, but excluded rows are shaded and marked
    Exclude_from_Excel = True.
    """
    global result_df, controller

    if result_df is None or result_df.empty:
        messagebox.showwarning("Exclude", "No results table is available yet.")
        return

    if controller is None or controller.tree is None:
        messagebox.showwarning("Exclude", "No table row is selected.")
        return

    selected = controller.tree.selection()

    if not selected:
        messagebox.showwarning("Exclude", "Please select a row in the table first.")
        return

    # Toggle all selected rows.
    for item in selected:
        vals = controller.tree.item(item, "values")

        if not vals:
            continue

        test_name = vals[0]
        mask = result_df["Test_Name"] == test_name

        if mask.any():
            current_value = bool(result_df.loc[mask, "Exclude_from_Excel"].iloc[0])
            new_value = not current_value
            result_df.loc[mask, "Exclude_from_Excel"] = new_value
            for idx in result_df.index[mask]:
                manual_exclusion_overrides[idx] = new_value

    update_table(result_df, controller)


def select_directory():
    d = filedialog.askdirectory()
    if d:
        entry_dir.delete(0, tk.END)
        entry_dir.insert(0, d)
        selected_dir_var.set(d)


def process_all():
    """Parse, analyze, plot, then populate table and wire controller."""
    global result_df, controller

    base_dir = selected_dir_var.get().strip()

    if not base_dir or not os.path.isdir(base_dir):
        messagebox.showwarning("Input", "Please select a valid directory.")
        return

    try:
        radius_nm = float(entry_radius.get().strip())
    except ValueError:
        messagebox.showwarning("Input", "Please enter a numeric radius (nm).")
        return

    # Validate the optional filter before running the expensive Hertz fits.
    try:
        get_r2_cutoff()
    except ValueError as exc:
        messagebox.showwarning("R² cutoff", str(exc))
        return

    items = [x for x in os.listdir(base_dir) if x not in ("Figure", "default", "Calib") and x.lower() != "calibration"]
    items = [x for x in items if os.path.isdir(os.path.join(base_dir, x))]

    names = []
    E_list = []
    plast_list = []
    plast_pct_list = []
    fit_r2_list = []
    fit_rmse_list = []
    fit_nrmse_list = []
    fit_start_idx_list = []
    fit_end_idx_list = []
    fit_direction_list = []

    for sub in items:
        recs = read_file(os.path.join(base_dir, sub))

        if not recs:
            continue

        for rec in recs:
            D = Data(rec["crosshead"], rec["load"], rec["times"], radius_nm, rec["filename"], base_dir)

            E_Pa, popt, fit_idx, fit_r2, fit_rmse_uN, fit_nrmse_pct = D.hertz_fit_linearized_loading_with_C(
                w_min=45,
                w_max_frac=1.0,
                r2_floor=0.0,
                force_low_frac=0.0,
                force_high_frac=1.0,
                x0_bound_span_factor=10.0,
                C_bound_span_factor=10.0,
            )

            # Poisson correction from your current code.
            E_corrected = E_Pa * (1 - 0.42 ** 2) if np.isfinite(E_Pa) else np.nan

            p1, p2, plast_nm = D.plasticity_points()
            max_disp_nm = float(np.nanmax(D.x_nm)) if D.x_nm.size else np.nan

            if np.isfinite(plast_nm) and np.isfinite(max_disp_nm) and max_disp_nm > 0:
                plast_pct = 100.0 * plast_nm / max_disp_nm
            else:
                plast_pct = np.nan

            # Plots.
            D.plot_time_vs_load()
            D.plot_time_vs_crosshead()
            D.plot_crosshead_vs_load_full_with_segment_fit(
                popt=popt,
                fit_idx=fit_idx,
                plast_pts=(p1, p2, plast_nm),
                fit_r2=fit_r2,
                fit_rmse_uN=fit_rmse_uN,
                fit_nrmse_pct=fit_nrmse_pct,
            )

            names.append(D.name)
            E_list.append(E_corrected if np.isfinite(E_corrected) else np.nan)
            plast_list.append(plast_nm if np.isfinite(plast_nm) else np.nan)
            plast_pct_list.append(plast_pct)
            fit_r2_list.append(fit_r2)
            fit_rmse_list.append(fit_rmse_uN)
            fit_nrmse_list.append(fit_nrmse_pct)
            fit_start_idx_list.append(fit_idx[0] if fit_idx is not None else np.nan)
            fit_end_idx_list.append(fit_idx[1] if fit_idx is not None else np.nan)
            fit_direction_list.append(popt[3] if popt is not None and len(popt) >= 4 else np.nan)

    result_df = pd.DataFrame({
        "Test_Name": names,
        "Exclude_from_Excel": [False] * len(names),
        "E_corrected_Pa": E_list,
        "Hertz_fit_R2": fit_r2_list,
        "Hertz_fit_RMSE_uN": fit_rmse_list,
        "Hertz_fit_NRMSE_percent": fit_nrmse_list,
        "Fit_Start_Index": fit_start_idx_list,
        "Fit_End_Index": fit_end_idx_list,
        "Fit_Direction": fit_direction_list,
        "Plasticity_nm": plast_list,
        "Plasticity_percent": plast_pct_list,
    })

    # Starting a new analysis clears overrides from the previous dataset.
    manual_exclusion_overrides.clear()
    apply_r2_cutoff(refresh_table=False)

    controller = Controller(data_root=base_dir)
    controller.test_names = names[:]

    controller.set_plot_dir("crosshead VS load")
    update_table(result_df, controller)


def update_table(df, controller):
    for w in right_frame.winfo_children():
        w.destroy()

    cols = list(df.columns)
    tree = ttk.Treeview(right_frame, columns=cols, show="headings", selectmode="extended")
    width = right_frame.winfo_width() or 800
    col_w = max(120, width // max(1, len(cols)))

    for c in cols:
        tree.heading(c, text=c)
        tree.column(c, width=col_w, anchor="center")

    if "Test_Name" in cols:
        tree.column("Test_Name", width=250, anchor="center")

    if "Exclude_from_Excel" in cols:
        tree.column("Exclude_from_Excel", width=150, anchor="center")

    tree.tag_configure("excluded", background="#dddddd", foreground="#777777")

    for _, row in df.iterrows():
        values = list(row)
        is_excluded = False

        if "Exclude_from_Excel" in df.columns:
            is_excluded = bool(row["Exclude_from_Excel"])

        if is_excluded:
            tree.insert("", "end", values=values, tags=("excluded",))
        else:
            tree.insert("", "end", values=values)

    yscroll = ttk.Scrollbar(right_frame, orient="vertical", command=tree.yview)
    xscroll = ttk.Scrollbar(right_frame, orient="horizontal", command=tree.xview)
    tree.configure(yscrollcommand=yscroll.set, xscrollcommand=xscroll.set)

    yscroll.pack(side="right", fill="y")
    xscroll.pack(side="bottom", fill="x")
    tree.pack(fill="both", expand=True)

    controller.tree = tree
    tree.bind("<Double-1>", controller.on_row_double_click)


def set_plot_dir_crosshead_load():
    controller.set_plot_dir("crosshead VS load")


def set_plot_dir_time_load():
    controller.set_plot_dir("time VS load")


def set_plot_dir_time_crosshead():
    controller.set_plot_dir("time VS crosshead")


# ------------------------------- Tk UI Setup --------------------------------

root = tk.Tk()
root.title("MTS Data Analyzer - Full Loading Hertz Fit")

try:
    root.state("zoomed")
except Exception:
    root.geometry("1600x900")

selected_dir_var = tk.StringVar(value="")

# Layout tuned: widen left plot area; push right panel further right.
top_frame = tk.Frame(root)
top_frame.place(x=900, y=0, width=700, height=200)

dir_frame = tk.Frame(top_frame)
dir_frame.pack(pady=5)

tk.Label(dir_frame, text="Select directory:").pack(side=tk.LEFT)
entry_dir = tk.Entry(dir_frame, width=50)
entry_dir.pack(side=tk.LEFT, padx=5)
tk.Button(dir_frame, text="Browse", command=select_directory).pack(side=tk.LEFT)

tk.Label(top_frame, text="Enter probe/sample radius (nm):").pack()
entry_radius = tk.Entry(top_frame)
entry_radius.pack()

# Optional R² export filter; blank retains every fit by default.
r2_frame = tk.Frame(top_frame)
r2_frame.pack(pady=(5, 0))
tk.Label(r2_frame, text="Auto-exclude fits with R² below:").pack(side=tk.LEFT, padx=(0, 5))
entry_r2_cutoff = tk.Entry(r2_frame, width=9)
entry_r2_cutoff.pack(side=tk.LEFT)
tk.Button(r2_frame, text="Apply R²", command=apply_r2_cutoff).pack(side=tk.LEFT, padx=5)
entry_r2_cutoff.bind("<Return>", lambda _event: apply_r2_cutoff())

btn_frame = tk.Frame(top_frame)
btn_frame.pack(pady=8)
tk.Button(btn_frame, text="Run", command=process_all).grid(row=0, column=0, padx=5)
tk.Button(btn_frame, text="Exclude/Include Selected", command=toggle_exclude_selected).grid(row=0, column=1, padx=5)
tk.Button(btn_frame, text="Export to Excel", command=export_to_excel).grid(row=0, column=2, padx=5)

# Left plot area.
left_frame = tk.Frame(root)
left_frame.place(x=0, y=200, width=900, height=650)

# Left-top plot selectors.
left_top = tk.Frame(root)
left_top.place(x=120, y=100, width=760, height=100)
tk.Button(left_top, text="Crosshead VS Load", width=20, command=set_plot_dir_crosshead_load).grid(row=0, column=0, pady=(3, 30), padx=6)
tk.Button(left_top, text="Time VS Load", width=20, command=set_plot_dir_time_load).grid(row=0, column=1, pady=(3, 30), padx=6)
tk.Button(left_top, text="Time VS Crosshead", width=20, command=set_plot_dir_time_crosshead).grid(row=0, column=2, pady=(3, 30), padx=6)
tk.Button(left_top, text="Previous", width=14, command=lambda: controller.prev_test()).grid(row=2, column=0, padx=6)
tk.Button(left_top, text="Next", width=14, command=lambda: controller.next_test()).grid(row=2, column=2, padx=6)

# Unified right panel.
right_frame = tk.Frame(root)
right_frame.place(x=980, y=200, width=600, height=650)

result_df = None
controller = Controller(data_root="")

root.mainloop()
