"""Compare frozen NetMHCpan scores with natural ovalbumin binding assays."""
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

from common import DATA, RESULTS, plot_style, save_summary

ASSAY = "dissociation constant KD (~IC50)|purified MHC/competitive/radioactivity|nM"


def prepare(assays, scores, sequence):
    """Require exact source matching and retain a single non-conflicting label per pair."""
    if set(assays.source) != {"UNIPROT:P01012"} or set(assays.assay) != {ASSAY}:
        raise ValueError("Unexpected source antigen or assay type")
    if not assays.peptide.map(lambda peptide: sequence.count(peptide) == 1).all():
        raise ValueError("The source sequence must contain each peptide exactly once")
    keys = ["allele", "peptide"]
    assays = assays.assign(label=assays.qualitative_measure.str.startswith("Positive").astype(int))
    groups = assays.groupby(keys).label.agg(["min", "max", "size"]).reset_index()
    consensus = groups[groups["min"].eq(groups["max"])].rename(columns={"min": "label"})
    audited = consensus.drop(columns="max").merge(scores, on=keys, how="left", validate="one_to_one")
    if not np.isfinite(scores[["ba_rank", "el_rank"]].to_numpy()).any():
        raise ValueError("No valid predictor ranks")
    return audited, int((groups["min"] != groups["max"]).sum())


def evaluate(table):
    rows = []
    for allele, frame in table.groupby("allele"):
        expected_length = 15 if allele.startswith("HLA") else 9
        complete = frame.dropna(subset=["ba_rank", "el_rank"]).copy()
        complete["length_deviation"] = (complete.peptide.str.len() - expected_length).abs()
        for name, column in (
            ("Binding affinity (BA)", "ba_rank"), ("Eluted ligand (EL)", "el_rank"),
            ("Length baseline", "length_deviation"),
        ):
            y = complete.label
            p = -complete[column]
            paired = int(len(y))
            supported = y.nunique() == 2
            rows.append({
                "allele": allele, "predictor": name, "observed_pairs": len(frame),
                "scored_pairs": paired, "positive": int(y.sum()), "negative": int(paired - y.sum()),
                "prevalence": float(y.mean()) if paired else np.nan,
                "roc_auc": roc_auc_score(y, p) if supported else np.nan,
                "average_precision": average_precision_score(y, p) if supported else np.nan,
            })
    return pd.DataFrame(rows)


def main():
    seq = "".join((DATA / "ova.fasta").read_text().splitlines()[1:])
    assays = pd.read_csv(DATA / "binding_assays.csv")
    scores = pd.read_csv(DATA / "external_scores.csv")
    table, conflicts = prepare(assays, scores, seq)
    metrics = evaluate(table)
    RESULTS.mkdir(exist_ok=True)
    metrics.to_csv(RESULTS / "metrics.csv", index=False, float_format="%.6f")

    sensitivity = []
    for row in table.dropna(subset=["ba_rank", "el_rank"]).itertuples():
        for method in ("ba", "el"):
            rank = getattr(row, method + "_rank")
            sensitivity.append({"allele": row.allele, "label": row.label, "method": method,
                                "rank_le_2": rank <= 2})
    counts = pd.DataFrame(sensitivity).groupby(["allele", "method", "label"]).rank_le_2.agg(
        ["sum", "count"]
    ).rename(columns={"sum": "rank_le_2", "count": "assayed"}).reset_index()
    counts.to_csv(RESULTS / "rank_threshold_diagnostics.csv", index=False)
    coverage = metrics.drop_duplicates("allele")[["allele", "observed_pairs", "scored_pairs",
                                                    "positive", "negative"]]
    coverage.to_csv(RESULTS / "evaluation_coverage.csv", index=False)

    plt = plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.3), constrained_layout=True)
    colors = ["#227f74", "#456fb0", "#a47a46"]
    for ax, measure, label in zip(axes, ("roc_auc", "average_precision"), ("AUROC", "Average precision")):
        panel = metrics.pivot(index="allele", columns="predictor", values=measure)
        panel.plot.bar(ax=ax, color=colors, legend=False)
        ax.set(xlabel="", ylabel=label, ylim=(0, 1.05), title=label + " on matched assay records")
        ax.tick_params(axis="x", labelrotation=18, labelsize=8)
        if measure == "roc_auc":
            ax.axhline(0.5, linestyle=":", color="#8d9298")
    axes[1].legend(frameon=False, fontsize=8, loc="upper left")
    fig.savefig(RESULTS / "performance.png")
    plt.close(fig)
    save_summary({
        "assay_rows": len(assays), "nonconflicting_pairs": len(table), "conflicting_pairs": conflicts,
        "pairs_scored_by_both_models": int(table[["ba_rank", "el_rank"]].notna().all(axis=1).sum()),
        "alleles_with_both_labels": sorted(metrics.loc[metrics.roc_auc.notna(), "allele"].unique().tolist()),
    })


if __name__ == "__main__":
    main()
