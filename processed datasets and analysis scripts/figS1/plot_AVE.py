from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ============================================================
# User settings
# ============================================================
DATA_DIR = Path("data")
OUTPUT_DIR = Path("figures")
OUTPUT_NAME = "regional_ssh_ke_salinity_temperature_2014_2018"

DATA_START_DATE = pd.Timestamp("2013-01-01")
DATA_END_DATE = pd.Timestamp("2018-12-31")
PLOT_START_DATE = pd.Timestamp("2014-01-01")
PLOT_END_DATE = pd.Timestamp("2018-12-31")

# "auto" supports both a 365-day LICOM calendar and Gregorian daily data.
DAILY_CALENDAR = "auto"       # "auto", "noleap", or "gregorian"

# SSH input is assumed to be in metres and is displayed in centimetres.
SSH_SCALE = 1
KE_SCALE = 1.0e17

# SSH anomaly reference. "jan2014" subtracts each experiment's January 2014
# monthly mean from both its daily and monthly series. "mean2014_2018" instead
# subtracts its five-year mean.
SSH_ANOMALY_REFERENCE = "jan2014"

EXPERIMENTS = {
    "eg": {"label": "Reg-EG", "color": "#1f77b4"},
    "jg": {"label": "Reg-JG", "color": "#d62728"},
    "jl": {"label": "Reg-JL", "color": "#2ca02c"},
}

# Typical journal double-column width is about 7.0-7.2 inches.
FIGSIZE = (7.2, 3.0)
DPI = 300
Y_PADDING_FRACTION = 0.075
THIN_LINEWIDTH = 0.55
THICK_LINEWIDTH = 1.65


PANEL_CONFIGS = [
    {
        "letter": "a",
        "filename": "average-ssh-{exp}.txt",
        "ylabel": "SSH (m)",
        "trend_unit": "cm yr$^{-1}$",
        "scale": SSH_SCALE,
        "ssh_anomaly": True,
    },
    {
        "letter": "b",
        "filename": "regional-integrated-total-ke-{exp}.txt",
        "ylabel": "Total KE ($10^{17}$ J)",
        "trend_unit": "$10^{17}$ J yr$^{-1}$",
        "scale": KE_SCALE,
        "ssh_anomaly": False,
    },
    {
        "letter": "c",
        "filename": "average-ss-{exp}.txt",
        "ylabel": "Salinity (psu)",
        "trend_unit": "psu yr$^{-1}$",
        "scale": 1.0,
        "ssh_anomaly": False,
    },
    {
        "letter": "d",
        "filename": "average-st-{exp}.txt",
        "ylabel": "Temperature (°C)",
        "trend_unit": "°C yr$^{-1}$",
        "scale": 1.0,
        "ssh_anomaly": False,
    },
]


plt.rcParams.update(
    {
        "font.size": 8,
        "axes.labelsize": 8.5,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


def load_daily_series(filepath):
    """Load one-column daily data and mask fill values."""
    if not filepath.is_file():
        raise FileNotFoundError(f"Input file not found: {filepath}")

    values = np.asarray(np.loadtxt(filepath), dtype=np.float64).reshape(-1)
    values[np.abs(values) >= 1.0e30] = np.nan
    return values


def build_daily_dates(number_of_values, filepath):
    """Build six-year dates and infer Gregorian or no-leap calendar."""
    gregorian_dates = pd.date_range(DATA_START_DATE, DATA_END_DATE, freq="D")
    noleap_dates = gregorian_dates[
        ~((gregorian_dates.month == 2) & (gregorian_dates.day == 29))
    ]

    if DAILY_CALENDAR == "auto":
        if number_of_values == len(noleap_dates):
            return noleap_dates
        if number_of_values == len(gregorian_dates):
            return gregorian_dates
        raise ValueError(
            f"{filepath} contains {number_of_values} daily values; expected "
            f"{len(noleap_dates)} (no-leap) or {len(gregorian_dates)} "
            "(Gregorian) for 2013-2018."
        )

    daily_dates = {
        "noleap": noleap_dates,
        "gregorian": gregorian_dates,
    }.get(DAILY_CALENDAR)

    if daily_dates is None:
        raise ValueError(
            "DAILY_CALENDAR must be 'auto', 'noleap', or 'gregorian'"
        )
    if number_of_values != len(daily_dates):
        raise ValueError(
            f"{filepath} contains {number_of_values} values, but "
            f"{DAILY_CALENDAR} requires {len(daily_dates)}."
        )
    return daily_dates


def monthly_mean(daily_dates, daily_values):
    """Calculate monthly arithmetic means of daily regional values."""
    return pd.Series(daily_values, index=daily_dates).resample("MS").mean()


def linear_trend_per_year(dates, values):
    """Return a linear trend in plotted units per year."""
    elapsed_years = (
        (dates - PLOT_START_DATE).days.to_numpy(dtype=np.float64) / 365.2425
    )
    valid = np.isfinite(values)
    if np.count_nonzero(valid) < 2:
        return np.nan
    return np.polyfit(elapsed_years[valid], values[valid], 1)[0]


def format_trend(value):
    if not np.isfinite(value):
        return "nan"
    magnitude = abs(value)
    if magnitude != 0.0 and (magnitude < 1.0e-3 or magnitude >= 1.0e3):
        return f"{value:+.2e}"
    return f"{value:+.3f}"


def ssh_reference(monthly_values):
    """Return the configured SSH reference in the plotted unit (cm)."""
    if SSH_ANOMALY_REFERENCE == "jan2014":
        reference = monthly_values.loc[pd.Timestamp("2014-01-01")]
    elif SSH_ANOMALY_REFERENCE == "mean2014_2018":
        reference = monthly_values.loc[
            PLOT_START_DATE:pd.Timestamp("2018-12-01")
        ].mean()
    else:
        raise ValueError(
            "SSH_ANOMALY_REFERENCE must be 'jan2014' or 'mean2014_2018'"
        )

    if not np.isfinite(reference):
        raise ValueError("The selected SSH anomaly reference is not finite")
    return float(reference)


def prepare_panel_data(config, plot_monthly_dates):
    """Load six years, form monthly means, and retain only 2014-2018."""
    panel_data = {}

    for experiment, style in EXPERIMENTS.items():
        filepath = DATA_DIR / config["filename"].format(exp=experiment)
        all_daily_values = load_daily_series(filepath) / config["scale"]
        all_daily_dates = build_daily_dates(len(all_daily_values), filepath)
        all_monthly_values = monthly_mean(all_daily_dates, all_daily_values)

        if len(all_monthly_values) != 72:
            raise ValueError(
                f"{filepath} generated {len(all_monthly_values)} months, "
                "but 72 months are required for 2013-2018."
            )

        if config["ssh_anomaly"]:
            reference = 0#ssh_reference(all_monthly_values)
            all_daily_values = all_daily_values - reference
            all_monthly_values = all_monthly_values - reference

        daily_mask = (
            (all_daily_dates >= PLOT_START_DATE)
            & (all_daily_dates <= PLOT_END_DATE)
        )
        daily_dates = all_daily_dates[daily_mask]
        daily_values = all_daily_values[daily_mask]
        plot_monthly_values = all_monthly_values.reindex(plot_monthly_dates)

        if plot_monthly_values.isna().all():
            raise ValueError(f"No valid 2014-2018 monthly values in {filepath}")

        monthly_values = plot_monthly_values.to_numpy(dtype=np.float64)
        panel_data[experiment] = {
            "style": style,
            "daily_dates": daily_dates,
            "daily": daily_values,
            "monthly": monthly_values,
            "trend": linear_trend_per_year(
                plot_monthly_dates,
                monthly_values,
            ),
        }

    return panel_data


def apply_padded_y_limits(ax, panel_data):
    """Cover every daily/monthly curve with 7.5% vertical padding."""
    arrays = []
    for result in panel_data.values():
        arrays.extend((result["daily"], result["monthly"]))

    finite_arrays = [array[np.isfinite(array)] for array in arrays]
    finite_arrays = [array for array in finite_arrays if array.size]
    if not finite_arrays:
        return

    combined = np.concatenate(finite_arrays)
    data_min = np.min(combined)
    data_max = np.max(combined)
    data_range = data_max - data_min

    if data_range == 0.0:
        padding = max(abs(data_min) * Y_PADDING_FRACTION, 1.0e-6)
    else:
        padding = data_range * Y_PADDING_FRACTION

    ax.set_ylim(data_min - padding, data_max + padding)


def draw_panel(ax, config, dates, panel_data, show_legend=False):
    """Draw one panel without an individual title."""
    for result in panel_data.values():
        style = result["style"]
        # ax.plot(
        #     result["daily_dates"],
        #     result["daily"],
        #     color=style["color"],
        #     linewidth=THIN_LINEWIDTH,
        #     alpha=0.25,
        #     zorder=2,
        # )
        ax.plot(
            dates,
            result["monthly"],
            color=style["color"],
            linewidth=THICK_LINEWIDTH,
            label=style["label"],
            zorder=3,
        )

    ax.text(
        0.015,
        0.975,
        f"({config['letter']})",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        fontweight="bold",
        color="black",
        zorder=6,
    )

    # Trends based on the 60 monthly means from 2014-2018.
    # for row, result in enumerate(panel_data.values()):
    #     style = result["style"]
    #     ax.text(
    #         0.985,
    #         0.965 - row * 0.09,
    #         f"{style['label']}: {format_trend(result['trend'])} "
    #         f"{config['trend_unit']}",
    #         transform=ax.transAxes,
    #         ha="right",
    #         va="top",
    #         fontsize=8,
    #         color=style["color"],
    #         bbox={
    #             "facecolor": "white",
    #             "edgecolor": "none",
    #             "alpha": 0.68,
    #             "pad": 0.8,
    #         },
    #         zorder=5,
    #     )

    ax.set_ylabel(config["ylabel"])
    ax.set_xlim(PLOT_START_DATE, PLOT_END_DATE)
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax.xaxis.set_minor_locator(mdates.MonthLocator(interval=3))
    ax.grid(axis="y", color="0.88", linewidth=0.5)
    ax.tick_params(direction="out", length=3)
    apply_padded_y_limits(ax, panel_data)

    if show_legend:
        ax.legend(
            loc="lower right",
            ncol=1,
            frameon=False,
            handlelength=2.2,
            borderaxespad=0.3,
        )


def main():
    plot_monthly_dates = pd.date_range(
        PLOT_START_DATE,
        "2018-12-01",
        freq="MS",
    )
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    fig, axs = plt.subplots(
        nrows=2,
        ncols=2,
        figsize=FIGSIZE,
        sharex=True,
    )

    for panel_index, (ax, config) in enumerate(zip(axs.flat, PANEL_CONFIGS)):
        panel_data = prepare_panel_data(config, plot_monthly_dates)
        draw_panel(
            ax,
            config,
            plot_monthly_dates,
            panel_data,
            show_legend=(panel_index == len(PANEL_CONFIGS) - 1),
        )

    # Explicitly suppress duplicated x-axis labels on the upper row.
    for ax in axs[0, :]:
        ax.tick_params(axis="x", which="both", labelbottom=False)

    for ax in axs[1, :]:
        ax.set_xlabel("Year")

    fig.subplots_adjust(
        left=0.105,
        right=0.99,
        bottom=0.105,
        top=0.985,
        wspace=0.27,
        hspace=0.16,
    )

    png_path = OUTPUT_DIR / f"{OUTPUT_NAME}.png"
    pdf_path = OUTPUT_DIR / f"{OUTPUT_NAME}.pdf"
    fig.savefig(png_path, dpi=DPI, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {png_path}")
    print(f"Saved: {pdf_path}")


if __name__ == "__main__":
    main()