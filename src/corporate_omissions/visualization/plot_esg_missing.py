"""Plotting helpers for ESG missing-data analysis (emissions decomposition)."""
from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D


class MissingESGData_Plots:
    """
    A class for plotting article graphs related to ESG data.
    """

    def __init__(self):
        pass

    # ------------------------------------------------------------------
    # Scope-1 disclosure composition (year / sector / size)
    # ------------------------------------------------------------------
    def graph_scope_1_disclosure_pct(self, raws, for_plot):
        """
        Generate stacked bar plots showing the proportion of Scope 1
        disclosures (Disclosed, Estimated, Missing) by year, sector,
        and firm size.

        Parameters
        ----------
        raws : object
            Must have a ``df_funds`` attribute containing the ESG data.
        for_plot : pd.DataFrame
            Processed data for plotting.
        """
        if not hasattr(raws, "df_funds"):
            raise AttributeError(
                "'raws' object must have a 'df_funds' attribute containing the ESG data."
            )

        # === 1. Plot by Year ===
        total_sample_year = raws.df_funds.groupby("year").gvkey.count()
        estimated_year = (
            raws.df_funds[~raws.df_funds.sc1_estimated.isna()]
            .groupby("year")
            .gvkey.count()
        )
        disclosed_year = (
            raws.df_funds[~raws.df_funds.sc1_disclosed.isna()]
            .groupby("year")
            .gvkey.count()
        )
        remaining_year = total_sample_year - estimated_year - disclosed_year

        df_stacked_year = (
            pd.DataFrame(
                {
                    "Disclosed": disclosed_year,
                    "Estimated": estimated_year,
                    "Remaining": remaining_year,
                }
            )
            .fillna(0)
            .reset_index()
        )

        # === 2. Plot by Sector ===
        disclosed = (
            for_plot.dropna(subset=["sc1_disclosed"]).groupby("gsector").gvkey.count()
            / for_plot.groupby("gsector").gvkey.count()
        )
        estimated = (
            for_plot.dropna(subset=["sc1_estimated"]).groupby("gsector").gvkey.count()
            / for_plot.groupby("gsector").gvkey.count()
        )
        remaining = 1 - disclosed - estimated

        df_stacked_sector = pd.DataFrame(
            {
                "gsector": disclosed.index,
                "Disclosed": disclosed,
                "Estimated": estimated,
                "Remaining": remaining,
            }
        ).fillna(0)

        df_stacked_sector[["Disclosed", "Estimated", "Remaining"]] *= 100
        df_stacked_sector = df_stacked_sector.round(2)

        # === 3. Plot by Firm Size Decile ===
        for_plot = for_plot.copy()
        for_plot["size_decile"] = pd.qcut(for_plot["at"], 10, labels=False)

        disclosed_size = (
            for_plot.dropna(subset=["sc1_disclosed"])
            .groupby("size_decile")
            .gvkey.count()
            / for_plot.groupby("size_decile").gvkey.count()
        )
        estimated_size = (
            for_plot.dropna(subset=["sc1_estimated"])
            .groupby("size_decile")
            .gvkey.count()
            / for_plot.groupby("size_decile").gvkey.count()
        )
        disclosed_size = disclosed_size.fillna(0)
        estimated_size = estimated_size.fillna(0)
        remaining_size = 1 - disclosed_size - estimated_size

        df_stacked_size = pd.DataFrame(
            {
                "Decile": disclosed_size.index,
                "Disclosed": disclosed_size,
                "Estimated": estimated_size,
                "Remaining": remaining_size,
            }
        ).fillna(0)

        df_stacked_size[["Disclosed", "Estimated", "Remaining"]] *= 100
        df_stacked_size = df_stacked_size.round(2)

        # === Plotting ===
        tab20_colors = plt.cm.tab20.colors[:2]
        remaining_color = "white"

        fig, axes = plt.subplots(3, 1, figsize=(12, 12), sharex=False)

        # Plot 1: By Year
        axes[0].bar(
            df_stacked_year["year"],
            df_stacked_year["Disclosed"],
            label="Disclosed",
            color=tab20_colors[0],
        )
        axes[0].bar(
            df_stacked_year["year"],
            df_stacked_year["Estimated"],
            bottom=df_stacked_year["Disclosed"],
            color=tab20_colors[1],
            label="Estimated",
        )
        axes[0].bar(
            df_stacked_year["year"],
            df_stacked_year["Remaining"],
            bottom=df_stacked_year["Disclosed"] + df_stacked_year["Estimated"],
            label="Missing",
            color=remaining_color,
            edgecolor="black",
            linestyle="--",
            linewidth=0.8,
        )
        axes[0].set_ylabel("Number of Firms")
        axes[0].set_title("Scope 1 Disclosure by Year")

        # Plot 2: By Sector
        axes[1].bar(
            df_stacked_sector["gsector"],
            df_stacked_sector["Disclosed"],
            label="Disclosed",
            color=tab20_colors[0],
        )
        axes[1].bar(
            df_stacked_sector["gsector"],
            df_stacked_sector["Estimated"],
            bottom=df_stacked_sector["Disclosed"],
            color=tab20_colors[1],
        )
        axes[1].bar(
            df_stacked_sector["gsector"],
            df_stacked_sector["Remaining"],
            bottom=df_stacked_sector["Disclosed"] + df_stacked_sector["Estimated"],
            color=remaining_color,
            edgecolor="black",
            linestyle="--",
            linewidth=0.8,
        )
        axes[1].set_ylabel("Proportion of Firms (%)")
        axes[1].set_title("Scope 1 Disclosure by Sector")
        axes[1].tick_params(axis="x", labelrotation=25)

        # Plot 3: By Firm Size
        axes[2].bar(
            df_stacked_size["Decile"],
            df_stacked_size["Disclosed"],
            label="Disclosed",
            color=tab20_colors[0],
        )
        axes[2].bar(
            df_stacked_size["Decile"],
            df_stacked_size["Estimated"],
            bottom=df_stacked_size["Disclosed"],
            color=tab20_colors[1],
        )
        axes[2].bar(
            df_stacked_size["Decile"],
            df_stacked_size["Remaining"],
            bottom=df_stacked_size["Disclosed"] + df_stacked_size["Estimated"],
            color=remaining_color,
            edgecolor="black",
            linestyle="--",
            linewidth=0.8,
        )
        axes[2].set_xlabel("Firm Size Decile (based on Market Cap)")
        axes[2].set_ylabel("Proportion of Firms (%)")
        axes[2].set_title("Scope 1 Disclosure by Firm Size")

        # Grid and Legend
        for ax in axes:
            ax.grid(axis="y", linestyle="--", alpha=0.5)
        axes[0].legend(
            title="Disclosure Type", bbox_to_anchor=(1, 1), loc="upper left"
        )
        plt.tight_layout()
        plt.show()

    # ------------------------------------------------------------------
    # Emissions decomposition (horizontal waterfall bars)
    # ------------------------------------------------------------------
    def plot_emissions(self, df, by_sector=False):
        """
        Plot emissions for a given year, either total or by sector.

        Parameters
        ----------
        df : pd.DataFrame
            DataFrame containing emissions data.
        by_sector : bool
            If True, plot emissions by sector. If False, plot total emissions.
        """
        colors = plt.cm.Paired.colors[:10]

        legend_handles = [
            mpatches.Patch(color=colors[2], label="Exact"),
            mpatches.Patch(color=colors[0], label="Estimated"),
            mpatches.Patch(color=colors[4], label="Missing"),
            Line2D([0], [0], color="black", linestyle="--", linewidth=1, label="Tilt"),
        ]

        if not by_sector:
            # === Total Emissions ===
            exact = df[~df.sc1_disclosed.isna()].sc1_adj.sum() / 1e6
            estimated = df[~df.sc1_estimated.isna()].sc1_estimated.sum() / 1e6
            estimated_tilts = (
                df[~df.sc1_estimated.isna()].sc1_adj.sum() / 1e6
                - df[~df.sc1_estimated.isna()].sc1_estimated.sum() / 1e6
            )
            imputed = df[df.sc1.isna()].sc1_disclosed_filled.sum() / 1e6
            imputed_tilts = (
                df[df.sc1.isna()].sc1_adj.sum() / 1e6
                - df[df.sc1.isna()].sc1_disclosed_filled.sum() / 1e6
            )

            exact_avg = df[~df.sc1_disclosed.isna()].sc1_adj.mean()
            estimated_avg = df[~df.sc1_estimated.isna()].sc1_estimated.mean()
            estimated_tilts_avg = (
                df[~df.sc1_estimated.isna()].sc1_adj.mean()
                - df[~df.sc1_estimated.isna()].sc1_estimated.mean()
            )
            imputed_avg = df[df.sc1.isna()].sc1_disclosed_filled.mean()
            imputed_tilts_avg = (
                df[df.sc1.isna()].sc1_adj.mean()
                - df[df.sc1.isna()].sc1_disclosed_filled.mean()
            )

            measures = [exact, estimated, estimated_tilts, imputed, imputed_tilts]
            measures = [np.round(m, 2) for m in measures]
            measures_avg = [
                exact_avg,
                estimated_avg,
                estimated_tilts_avg,
                imputed_avg,
                imputed_tilts_avg,
            ]
            measures_avg = [np.round(m, 6) for m in measures_avg]

            lefts = [sum(measures[:i]) for i in range(len(measures))]
            lefts_avg = [sum(measures_avg[:i]) for i in range(len(measures_avg))]
            fill_colors = [colors[2], colors[0], colors[0], colors[4], colors[4]]
            linestyles = ["solid", "solid", "dashed", "solid", "dashed"]

            # --- Total emissions bar ---
            fig, ax = plt.subplots(figsize=(20, 0.7))
            for i in range(len(measures)):
                bar = ax.barh(
                    0,
                    measures[i],
                    left=lefts[i],
                    color=fill_colors[i],
                    edgecolor="black",
                    linewidth=0.8,
                )
                if linestyles[i] == "dashed":
                    for patch in bar:
                        patch.set_linestyle("--")
                center_x = lefts[i] + measures[i] / 2
                print(measures[i])
                if i == 2:  # estimated tilt
                    text_x = lefts[i] + measures[i] - 1
                    text_y = 1.05
                    ax.annotate(
                        f"{measures[i]:,.0f}",
                        xy=(center_x, 0),
                        xytext=(text_x, text_y),
                        ha="right",
                        va="center",
                        fontsize=11,
                        arrowprops=dict(arrowstyle="-", color="black"),
                    )
                elif i == 1:
                    text_x = lefts[i] + measures[i] - 1
                    text_y = -1.05
                    ax.annotate(
                        f"{measures[i]:,.0f}",
                        xy=(center_x, 0),
                        xytext=(text_x, text_y),
                        ha="right",
                        va="center",
                        fontsize=11,
                        arrowprops=dict(arrowstyle="-", color="black"),
                    )
                else:
                    ax.text(
                        center_x,
                        0,
                        f"{measures[i]:,.0f}",
                        ha="center",
                        va="center",
                        fontsize=11,
                    )

            ax.axis("off")
            ax.legend(
                handles=legend_handles,
                loc="upper center",
                bbox_to_anchor=(0.5, -0.2),
                ncol=4,
                frameon=False,
            )
            plt.title("Total Emissions (Million Metric Tons CO2)", fontsize=12)
            plt.tight_layout()
            plt.show()

            total_emissions = sum(measures)
            print(
                f"Total Emissions (GT CO2): {np.round(total_emissions / 1000, 5)}"
            )
            print(
                f"Difference vs EPA: "
                f"{np.round((np.round(total_emissions, 5) - 6343.2) / 6343.2 * 100, 2)}%"
            )

            # --- Average emissions bar ---
            fig, ax = plt.subplots(figsize=(20, 0.7))
            for i in range(len(measures_avg)):
                bar = ax.barh(
                    0,
                    measures_avg[i],
                    left=lefts_avg[i],
                    color=fill_colors[i],
                    edgecolor="black",
                    linewidth=0.8,
                )
                if linestyles[i] == "dashed":
                    for patch in bar:
                        patch.set_linestyle("--")
                center_x = lefts_avg[i] + measures_avg[i] / 2
                print(measures_avg[i])
                if i == 2:  # estimated tilt
                    text_x = lefts_avg[i] + measures_avg[i] + 0.05
                    text_y = -1.55
                    ax.annotate(
                        f"{measures_avg[i]:,.0f}",
                        xy=(center_x, 0),
                        xytext=(text_x, text_y),
                        ha="right",
                        va="center",
                        fontsize=11,
                        arrowprops=dict(arrowstyle="-", color="black"),
                    )
                elif i == 1:
                    text_x = lefts_avg[i] + measures_avg[i] - 0.05
                    text_y = -1.05
                    ax.annotate(
                        f"{measures_avg[i]:,.0f}",
                        xy=(center_x, 0),
                        xytext=(text_x, text_y),
                        ha="right",
                        va="center",
                        fontsize=11,
                        arrowprops=dict(arrowstyle="-", color="black"),
                    )
                else:
                    ax.text(
                        center_x,
                        0,
                        f"{measures_avg[i]:,.0f}",
                        ha="center",
                        va="center",
                        fontsize=11,
                    )

            ax.axis("off")
            ax.legend(
                handles=legend_handles,
                loc="upper center",
                bbox_to_anchor=(0.5, -1.5),
                ncol=4,
                frameon=False,
            )
            plt.title("Average Emissions (Metric Tons CO2)", fontsize=12)
            plt.tight_layout()
            plt.show()

        else:
            # === Emissions by Sector ===
            df = df.copy()
            df["Exact"] = df[~df.sc1_disclosed.isna()].sc1_disclosed / 1e6
            df["Estimated"] = df[~df.sc1_estimated.isna()].sc1_adj / 1e6
            df["Estimated_Tilts"] = (
                df[~df.sc1_estimated.isna()].sc1_adj / 1e6
                - df[~df.sc1_estimated.isna()].sc1_estimated / 1e6
            )
            df["Imputed"] = df[df.sc1.isna()].sc1_disclosed_filled / 1e6
            df["Imputed_Tilts"] = (
                df[df.sc1.isna()].sc1_adj / 1e6
                - df[df.sc1.isna()].sc1_disclosed_filled / 1e6
            )

            sectors = df["gsector"].unique()
            sectors_data = {sector: [] for sector in sectors}

            for sector in sectors:
                sector_data = df[df["gsector"] == sector]
                sector_measures = [
                    sector_data["Exact"].sum(),
                    sector_data["Estimated"].sum(),
                ]
                if sector_data["Estimated_Tilts"].sum() > 0:
                    sector_measures.append(sector_data["Estimated_Tilts"].sum())
                if sector_data["Imputed_Tilts"].sum() > 0:
                    sector_measures.append(sector_data["Imputed_Tilts"].sum())
                sector_measures.append(sector_data["Imputed"].sum())
                sectors_data[sector] = sector_measures

            # Normalize
            for sector in sectors_data:
                total = sum(sectors_data[sector])
                sectors_data[sector] = (
                    [m / total for m in sectors_data[sector]]
                    if total != 0
                    else [0 for m in sectors_data[sector]]
                )

            fig, ax = plt.subplots(figsize=(20, 10))
            lefts_arr = np.zeros(len(sectors))
            fill_colors = [colors[2], colors[0], colors[0], colors[4], colors[4]]
            linestyles = ["solid", "solid", "dashed", "solid", "dashed"]

            for i, (sector, sector_measures) in enumerate(sectors_data.items()):
                for j in range(len(sector_measures)):
                    bar = ax.barh(
                        y=i,
                        width=sector_measures[j],
                        left=lefts_arr[i],
                        color=fill_colors[j],
                        edgecolor="black",
                        linewidth=0.8,
                    )
                    if linestyles[j] == "dashed":
                        for patch in bar:
                            patch.set_linestyle("--")
                    center_x = lefts_arr[i] + sector_measures[j] / 2
                    ax.text(
                        center_x,
                        i,
                        f"{sector_measures[j]:.2f}",
                        va="center",
                        ha="center",
                        fontsize=11,
                        color="black",
                    )
                    lefts_arr[i] += sector_measures[j]

            ax.axis("off")
            ax.legend(
                handles=legend_handles,
                loc="upper center",
                bbox_to_anchor=(0.5, 0),
                ncol=4,
                frameon=False,
            )
            ax.set_yticks(np.arange(len(sectors)))
            ax.set_yticklabels(sectors)

            for i, sector in enumerate(sectors):
                ax.text(
                    -0.05,
                    i,
                    sector,
                    va="center",
                    ha="right",
                    fontsize=10,
                    color="black",
                )

            plt.title("Standardized Emissions by Sector (Million Metric Tons CO2)")
            plt.tight_layout()
            plt.show()

    # ------------------------------------------------------------------
    # Emissions bars (observed vs tilted, total + average)
    # ------------------------------------------------------------------
    def plot_emissions_bars(self, df, by_sector=False):
        """
        Plot emissions for a given year, either total or by sector.

        Parameters
        ----------
        df : pd.DataFrame
            DataFrame containing emissions data.
        by_sector : bool
            If True, plot emissions by sector. If False, plot total emissions.
        """
        colors = plt.cm.Paired.colors[:10]

        legend_handles = [
            mpatches.Patch(color=colors[2], label="Exact"),
            mpatches.Patch(color=colors[0], label="Estimated"),
            mpatches.Patch(color=colors[4], label="Missing"),
            Line2D([0], [0], color="black", linestyle="--", linewidth=1, label="Tilt"),
        ]

        if not by_sector:
            # === Total Emissions ===
            exact = df[~df.sc1_disclosed.isna()].sc1_adj.sum() / 1e6
            estimated = df[~df.sc1_estimated.isna()].sc1_estimated.sum() / 1e6
            estimated_tilt = (
                df[~df.sc1_estimated.isna()].sc1_adj.sum() / 1e6 - estimated
            )
            missing = df[df.sc1.isna()].sc1_disclosed_filled.sum() / 1e6
            missing_tilt = df[df.sc1.isna()].sc1_adj.sum() / 1e6 - missing

            bar1_values = [exact, estimated, missing]
            bar1_colors = [colors[2], colors[0], colors[4]]

            bar2_values = [exact, estimated, estimated_tilt, missing, missing_tilt]
            bar2_colors = [colors[2], colors[0], colors[0], colors[4], colors[4]]
            bar2_linewidths = [0.8, 0.8, 1.5, 0.8, 1.5]

            fig, ax = plt.subplots(figsize=(4, 6))

            # Plot first bar (non-tilted)
            left = 0
            for value, color in zip(bar1_values, bar1_colors):
                ax.bar(
                    0.5,
                    value,
                    bottom=left,
                    width=0.15,
                    color=color,
                    edgecolor="black",
                    linewidth=0.8,
                )
                left += value

            # Plot second bar (tilted)
            left = 0
            for value, color, lw in zip(bar2_values, bar2_colors, bar2_linewidths):
                ax.bar(
                    0.8,
                    value,
                    bottom=left,
                    width=0.15,
                    color=color,
                    edgecolor="black",
                    linewidth=lw,
                )
                left += value

            # Add value labels on bars
            for bar_values, x in zip([bar1_values, bar2_values], [0.5, 0.8]):
                y_bottom = 0
                for v in bar_values:
                    if v > 0.1:
                        if v < 110:
                            ax.plot(
                                [x, x + 0.1],
                                [y_bottom + v / 2, y_bottom],
                                color="black",
                                lw=0.8,
                            )
                            ax.text(
                                x + 0.1,
                                y_bottom,
                                f"{v:,.0f}",
                                va="center",
                                ha="left",
                                fontsize=9,
                            )
                        elif v < 200:
                            bar_top = y_bottom + v
                            ax.plot(
                                [x, x + 0.1],
                                [y_bottom + v / 2, bar_top],
                                color="black",
                                lw=0.8,
                            )
                            ax.text(
                                x + 0.1,
                                bar_top,
                                f"{v:,.0f}",
                                va="center",
                                ha="left",
                                fontsize=9,
                            )
                        else:
                            ax.text(
                                x,
                                y_bottom + v / 2,
                                f"{v:,.0f}",
                                ha="center",
                                va="center",
                                fontsize=9,
                            )
                    y_bottom += v

            ax.axhline(
                y=6343.2, color="red", linestyle="--", linewidth=1, label="EPA Estimate"
            )

            legend_handles_total = [
                mpatches.Patch(color=colors[2], label="Exact"),
                Line2D(
                    [0], [0], color="black", linestyle="-", linewidth=1.5, label="Tilt"
                ),
                mpatches.Patch(color=colors[0], label="Estimated"),
                Line2D(
                    [0],
                    [0],
                    color="red",
                    linestyle="--",
                    linewidth=1,
                    label="EPA Estimate",
                ),
                mpatches.Patch(color=colors[4], label="Missing"),
            ]

            ax.legend(
                handles=legend_handles_total,
                loc="upper center",
                bbox_to_anchor=(0.5, -0.1),
                ncol=3,
                frameon=False,
            )
            ax.set_xticks([0.5, 0.8])
            ax.set_xticklabels(["Observed \n Distribution", "Tilted"])
            ax.set_xlim([0.3, 1])
            ax.set_ylabel("Emissions (Million Metric Tons CO\u2082)")
            ax.grid(alpha=0.3, linestyle="--")
            plt.title("Total Emissions with and without Tilts")
            plt.tight_layout()
            plt.show()

            total_emissions = sum(bar2_values)
            print(
                f"Total Emissions (GT CO2): {np.round(total_emissions / 1000, 5)}"
            )
            print(
                f"Difference vs EPA: "
                f"{np.round((total_emissions - 6343.2) / 6343.2 * 100, 2)}%"
            )

            # === Average Emissions ===
            exact_avg = df[~df.sc1_disclosed.isna()].sc1_adj.mean()
            estimated_avg = df[~df.sc1_estimated.isna()].sc1_estimated.mean()
            estimated_tilts_avg = (
                df[~df.sc1_estimated.isna()].sc1_adj.mean()
                - df[~df.sc1_estimated.isna()].sc1_estimated.mean()
            )
            missing_avg = df[df.sc1.isna()].sc1_disclosed_filled.mean()
            missing_tilts_avg = (
                df[df.sc1.isna()].sc1_adj.mean()
                - df[df.sc1.isna()].sc1_disclosed_filled.mean()
            )

            bar1_values = [exact_avg, estimated_avg, missing_avg]
            bar1_colors = [colors[2], colors[0], colors[4]]

            bar2_values = [
                exact_avg,
                estimated_avg,
                estimated_tilts_avg,
                missing_avg,
                missing_tilts_avg,
            ]
            bar2_colors = [colors[2], colors[0], colors[0], colors[4], colors[4]]
            bar2_linewidths = [0.8, 0.8, 1.5, 0.8, 1.5]

            fig, ax = plt.subplots(figsize=(4, 6))

            # Plot first bar (non-tilted)
            left = 0
            for value, color in zip(bar1_values, bar1_colors):
                ax.bar(
                    0.5,
                    value,
                    bottom=left,
                    width=0.15,
                    color=color,
                    edgecolor="black",
                    linewidth=0.8,
                )
                left += value

            # Plot second bar (tilted)
            left = 0
            for value, color, lw in zip(bar2_values, bar2_colors, bar2_linewidths):
                ax.bar(
                    0.8,
                    value,
                    bottom=left,
                    width=0.15,
                    color=color,
                    edgecolor="black",
                    linewidth=lw,
                )
                left += value

            # Add value labels on bars
            for bar_values, x in zip([bar1_values, bar2_values], [0.5, 0.8]):
                y_bottom = 0
                for v in bar_values:
                    if v > 0.1:
                        if v < 55000:
                            ax.plot(
                                [x, x + 0.1],
                                [y_bottom + v / 2, y_bottom],
                                color="black",
                                lw=0.8,
                            )
                            ax.text(
                                x + 0.1,
                                y_bottom,
                                f"{v:,.0f}",
                                va="center",
                                ha="left",
                                fontsize=9,
                            )
                        elif v < 60000:
                            bar_top = y_bottom + v
                            ax.plot(
                                [x, x + 0.1],
                                [y_bottom + v / 2, bar_top],
                                color="black",
                                lw=0.8,
                            )
                            ax.text(
                                x + 0.1,
                                bar_top,
                                f"{v:,.0f}",
                                va="center",
                                ha="left",
                                fontsize=9,
                            )
                        else:
                            ax.text(
                                x,
                                y_bottom + v / 2,
                                f"{v:,.0f}",
                                ha="center",
                                va="center",
                                fontsize=9,
                            )
                    y_bottom += v

            legend_handles_avg = [
                mpatches.Patch(color=colors[2], label="Exact"),
                mpatches.Patch(color=colors[0], label="Estimated"),
                mpatches.Patch(color=colors[4], label="Missing"),
                Line2D(
                    [0], [0], color="black", linestyle="-", linewidth=1.5, label="Tilt"
                ),
            ]

            ax.legend(
                handles=legend_handles_avg,
                loc="upper center",
                bbox_to_anchor=(0.5, -0.1),
                ncol=4,
                frameon=False,
            )
            ax.set_xticks([0.5, 0.8])
            ax.set_xticklabels(["Observed \n Distribution", "Tilted"])
            ax.set_xlim([0.3, 1])
            ax.set_ylabel("Emissions (Metric Tons CO\u2082)")
            ax.grid(alpha=0.3, linestyle="--")
            plt.title("Average Emissions with and without Tilts", pad=20)
            plt.tight_layout()
            plt.show()

        else:
            # === Emissions by Sector ===
            df = df.copy()
            df["Exact"] = df[~df.sc1_disclosed.isna()].sc1_disclosed / 1e6
            df["Estimated"] = df[~df.sc1_estimated.isna()].sc1_adj / 1e6
            df["Estimated_Tilts"] = (
                df[~df.sc1_estimated.isna()].sc1_adj / 1e6
                - df[~df.sc1_estimated.isna()].sc1_estimated / 1e6
            )
            df["Imputed"] = df[df.sc1.isna()].sc1_disclosed_filled / 1e6
            df["Imputed_Tilts"] = (
                df[df.sc1.isna()].sc1_adj / 1e6
                - df[df.sc1.isna()].sc1_disclosed_filled / 1e6
            )

            sectors = df["gsector"].unique()
            sectors_data = {sector: [] for sector in sectors}

            for sector in sectors:
                sector_data = df[df["gsector"] == sector]
                sector_measures = [
                    sector_data["Exact"].sum(),
                    sector_data["Estimated"].sum(),
                ]
                if sector_data["Estimated_Tilts"].sum() > 0:
                    sector_measures.append(sector_data["Estimated_Tilts"].sum())
                if sector_data["Imputed_Tilts"].sum() > 0:
                    sector_measures.append(sector_data["Imputed_Tilts"].sum())
                sector_measures.append(sector_data["Imputed"].sum())
                sectors_data[sector] = sector_measures

            # Normalize
            for sector in sectors_data:
                total = sum(sectors_data[sector])
                sectors_data[sector] = (
                    [m / total for m in sectors_data[sector]]
                    if total != 0
                    else [0 for m in sectors_data[sector]]
                )

            fig, ax = plt.subplots(figsize=(20, 10))
            bottoms = np.zeros(len(sectors))
            fill_colors = [colors[2], colors[0], colors[0], colors[4], colors[4]]
            linestyle_defs = ["solid", "solid", "dashed", "solid", "dashed"]

            for j in range(
                max(len(measures) for measures in sectors_data.values())
            ):
                heights = []
                for sector in sectors:
                    sector_measures = sectors_data[sector]
                    if j < len(sector_measures):
                        heights.append(sector_measures[j])
                    else:
                        heights.append(0)

                bar = ax.bar(
                    x=np.arange(len(sectors)),
                    height=heights,
                    bottom=bottoms,
                    color=fill_colors[j],
                    edgecolor="black",
                    linewidth=0.8,
                )
                if linestyle_defs[j] == "dashed":
                    for patch in bar:
                        patch.set_linestyle("-")
                        patch.set_linewidth(1.5)

                for i in range(len(sectors)):
                    if heights[i] > 0:
                        ax.text(
                            i,
                            bottoms[i] + heights[i] / 2,
                            f"{heights[i]:.2f}",
                            ha="center",
                            va="center",
                            fontsize=11,
                            color="black",
                        )
                    bottoms[i] += heights[i]

            ax.set_xticks(np.arange(len(sectors)))
            ax.set_xticklabels(sectors, rotation=45, ha="right")
            ax.set_ylim([0, 1.1])
            ax.set_ylabel("Standardized Emissions Share")
            ax.legend(
                handles=legend_handles,
                loc="upper center",
                bbox_to_anchor=(0.5, -0.1),
                ncol=4,
                frameon=False,
            )
            ax.grid(axis="y", alpha=0.3, linestyle="--")
            plt.title(
                "Standardized Emissions by Sector (Million Metric Tons CO\u2082)"
            )
            plt.tight_layout()
            plt.show()