import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


# to run use followiong format or run from bash cript with venv   python pca_rscu.py --config config.json --output pca_output.csv

def load_config(config_path):
    with open(config_path, "r", encoding="utf-8") as f:
        return json.load(f)


def detect_sample_id_column(df):
    candidates = ["id", "sample", "sample_id", "name", "ID", "Sample"]
    for col in candidates:
        if col in df.columns:
            return col
    return df.columns[0]


def validate_numeric(df, source_name):
    out = df.copy()

    for col in out.columns:
        out[col] = pd.to_numeric(out[col], errors="coerce")

    if out.isna().any().any():
        bad_cols = out.columns[out.isna().any()].tolist()
        raise ValueError(
            f"{source_name}: missing or non-numeric values found in columns: {bad_cols}"
        )

    return out


def parse_metadata(sample_name):
    parts = str(sample_name).split("_")
    strain = parts[0] if len(parts) >= 1 else "unknown"
    replicate = parts[1] if len(parts) >= 2 else "unknown"
    return strain, replicate


def load_wide_orientation(df, label, codon_index, drop_codons):
    """
    Format:
        rows = samples
        columns = codons

    Example:
        id,GCT,GCC,GCA,...
        68u201_rep1,1.1,0.9,1.0,...
        G7R_rep1,0.8,1.2,1.1,...
    """
    sample_id_col = detect_sample_id_column(df)

    df.columns = [str(c).strip() for c in df.columns]
    df = df.rename(
        columns={
            col: str(col).strip().upper()
            for col in df.columns
            if col != sample_id_col
        }
    )

    missing = [c for c in codon_index if c not in df.columns]
    if missing:
        raise ValueError(f"{label}: missing codon columns: {missing}")

    keep_codons = [c for c in codon_index if c not in drop_codons]

    metadata = pd.DataFrame()
    metadata["sample"] = df[sample_id_col].astype(str).str.strip()
    metadata["dataset"] = label
    metadata[["strain", "replicate"]] = metadata["sample"].apply(
        lambda s: pd.Series(parse_metadata(s))
    )

    codon_data = df[keep_codons].copy()
    codon_data = validate_numeric(codon_data, label)

    return codon_data, metadata


def load_long_orientation(df, label, codon_index, drop_codons):
    """
    Format:
        rows = codons
        columns = samples

    Example:
        codon,68u201_rep1,68u201_rep2,G7R_rep1
        GCT,1.1,1.0,0.8
        GCC,0.9,1.1,1.2
        GCA,1.0,0.9,1.1
    """
    codon_col_candidates = ["codon", "Codon", "CODON"]
    codon_col = None

    for col in codon_col_candidates:
        if col in df.columns:
            codon_col = col
            break

    if codon_col is None:
        codon_col = df.columns[0]

    df[codon_col] = df[codon_col].astype(str).str.strip().str.upper()
    df = df.set_index(codon_col)

    duplicates = df.index[df.index.duplicated()].unique().tolist()
    if duplicates:
        raise ValueError(f"{label}: duplicate codons found: {duplicates}")

    missing = [c for c in codon_index if c not in df.index]
    if missing:
        raise ValueError(f"{label}: missing codon rows: {missing}")

    keep_codons = [c for c in codon_index if c not in drop_codons]

    df = df.loc[keep_codons]

    # transpose so rows = samples and columns = codons
    codon_data = df.T.copy()
    codon_data = validate_numeric(codon_data, label)

    metadata = pd.DataFrame()
    metadata["sample"] = codon_data.index.astype(str)
    metadata["dataset"] = label
    metadata[["strain", "replicate"]] = metadata["sample"].apply(
        lambda s: pd.Series(parse_metadata(s))
    )

    codon_data = codon_data.reset_index(drop=True)

    return codon_data, metadata


def detect_orientation(df, codon_index):
    """
    Returns:
        "wide" = samples are rows, codons are columns
        "long" = codons are rows, samples are columns
    """
    normalized_columns = {str(c).strip().upper() for c in df.columns}
    codons_in_columns = sum(c in normalized_columns for c in codon_index)

    first_col = df.columns[0]
    first_col_values = set(df[first_col].astype(str).str.strip().str.upper())
    codons_in_first_column = sum(c in first_col_values for c in codon_index)

    if codons_in_columns >= 10:
        return "wide"

    if codons_in_first_column >= 10:
        return "long"

    raise ValueError(
        "Could not detect file orientation. Expected either codons as columns "
        "or codons in the first column."
    )


def load_codon_file(file_cfg, codon_index, drop_codons):
    label = file_cfg["label"]
    path = Path(file_cfg["path"])

    if not path.exists():
        raise FileNotFoundError(f"{label}: file not found: {path}")

    df = pd.read_csv(path)
    if df.empty:
        raise ValueError(f"{label}: input file is empty")

    codon_index = [c.upper() for c in codon_index]
    drop_codons = [c.upper() for c in drop_codons]

    orientation = file_cfg.get("orientation", "auto")

    if orientation == "auto":
        orientation = detect_orientation(df, codon_index)

    print(f"{label}: detected orientation = {orientation}")

    if orientation == "wide":
        return load_wide_orientation(df, label, codon_index, drop_codons)

    if orientation == "long":
        return load_long_orientation(df, label, codon_index, drop_codons)

    raise ValueError(
        f"{label}: orientation must be 'auto', 'wide', or 'long'. Got: {orientation}"
    )


def run_pca(codon_df, metadata_df):
    if codon_df.shape[0] < 2:
        raise ValueError("Need at least 2 samples to run PCA.")

    if codon_df.shape[1] < 2:
        raise ValueError("Need at least 2 codon features to run PCA.")

    x = StandardScaler().fit_transform(codon_df)

    pca = PCA(n_components=2)
    components = pca.fit_transform(x)

    pca_df = pd.DataFrame(
        components,
        columns=["principal component 1", "principal component 2"],
    )

    pca_df = pd.concat([metadata_df.reset_index(drop=True), pca_df], axis=1)

    return pca_df, pca


def plot_pca(pca_df, plot_cfg, pca_model):
    figsize = tuple(plot_cfg.get("figsize", [10, 10]))

    g = sns.relplot(
        data=pca_df,
        x="principal component 1",
        y="principal component 2",
        hue="strain",
        style="dataset",
        height=figsize[1],
        aspect=figsize[0] / figsize[1],
    )

    g.set_axis_labels(
        f"{plot_cfg.get('x_label', 'Principal Component - 1')} "
        f"({pca_model.explained_variance_ratio_[0]:.1%})",
        f"{plot_cfg.get('y_label', 'Principal Component - 2')} "
        f"({pca_model.explained_variance_ratio_[1]:.1%})",
    )

    g.fig.suptitle(plot_cfg.get("title", "PCA of RSCU data"), y=1.02)

    plt.show()


def main():
    parser = argparse.ArgumentParser(description="Run PCA on RSCU codon data")
    parser.add_argument("-c", "--config", default="config.json")
    parser.add_argument("-o", "--output", default=None)
    args = parser.parse_args()

    config = load_config(args.config)

    codon_tables = []
    metadata_tables = []

    for file_cfg in config["files"]:
        codon_data, metadata = load_codon_file(
            file_cfg,
            config["codons"]["index"],
            config["codons"]["drop"],
        )

        codon_tables.append(codon_data)
        metadata_tables.append(metadata)

    combined_codons = pd.concat(codon_tables, axis=0, ignore_index=True)
    combined_metadata = pd.concat(metadata_tables, axis=0, ignore_index=True)

    print("Combined codon matrix shape:", combined_codons.shape)
    print("Combined metadata shape:", combined_metadata.shape)

    pca_df, pca_model = run_pca(combined_codons, combined_metadata)

    print(pca_df)

    if args.output:
        pca_df.to_csv(args.output, index=False)
        print(f"Saved PCA output to {args.output}")

    plot_pca(pca_df, config["plot"], pca_model)


if __name__ == "__main__":
    main()
