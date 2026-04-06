import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def load_config(config_path: str) -> dict:
    with open(config_path, "r", encoding="utf-8") as f:
        return json.load(f)


def load_and_prepare_dataframe(file_cfg: dict, codon_index: list[str], drop_codons: list[str]) -> pd.DataFrame:
    df = pd.read_csv(file_cfg["path"])
    df.columns = [f'{file_cfg["label"]}_{col}' for col in df.columns]

    df = df.set_index(pd.Index(codon_index))
    df = df.drop(drop_codons, errors="ignore")
    df = df.iloc[:, 1:]   # keep behavior from your original script
    df = df.T

    return df


def parse_metadata(sample_name: str) -> tuple[str, str]:
    # expected format like: d1_68u201_...
    parts = sample_name.split("_")
    cell_line = parts[0] if len(parts) > 0 else ""
    strain = parts[1] if len(parts) > 1 else ""
    return cell_line, strain


def run_pca(df: pd.DataFrame) -> pd.DataFrame:
    x = StandardScaler().fit_transform(df)
    components = PCA(n_components=2).fit_transform(x)

    result = pd.DataFrame(
        components,
        columns=["principal component 1", "principal component 2"],
        index=df.index,
    )

    result["sample"] = result.index
    result[["Cell line", "strain"]] = result["sample"].to_series().apply(
        lambda s: pd.Series(parse_metadata(s))
    )

    return result


def plot_pca(pca_df: pd.DataFrame, plot_cfg: dict) -> None:
    figsize = tuple(plot_cfg.get("figsize", [10, 10]))

    g = sns.FacetGrid(pca_df, hue="strain", height=figsize[1])
    g.map(sns.scatterplot, "principal component 1", "principal component 2", alpha=0.7)
    g.add_legend()

    g.set_axis_labels(
        plot_cfg.get("x_label", "Principal Component - 1"),
        plot_cfg.get("y_label", "Principal Component - 2"),
    )
    g.fig.suptitle(plot_cfg.get("title", "PCA"), y=1.02)
    g.fig.set_size_inches(*figsize)

    plt.show()


def main():
    parser = argparse.ArgumentParser(description="Run PCA on RSCU data")
    parser.add_argument(
        "-c", "--config",
        default="config.json",
        help="Path to config JSON file"
    )
    args = parser.parse_args()

    config = load_config(args.config)

    codon_index = config["codons"]["index"]
    drop_codons = config["codons"]["drop"]

    dfs = [
        load_and_prepare_dataframe(file_cfg, codon_index, drop_codons)
        for file_cfg in config["files"]
    ]

    combined_df = pd.concat(dfs)
    print("Combined input:")
    print(combined_df)

    pca_df = run_pca(combined_df)

    print("\nPCA result:")
    print(pca_df)

    plot_pca(pca_df, config["plot"])


if __name__ == "__main__":
    main()
