import pandas as pd
import argparse
import logging
import os

def main() -> int:
    parser = argparse.ArgumentParser(description="Decompose syngraph table into per-node assignment TSVs")
    parser.add_argument("--syngtab", required=True, help="Path to the syngraph table (TSV)")
    parser.add_argument("--out-dir", required=True, help="Directory to write per-node files")
    args = parser.parse_args()

    # args.syngtab = 'data/hymenoptera/syngraph_run/hymenoptera.syngraph_tabulate.table.tsv'

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    if not os.path.exists(args.syngtab):
        logging.error("syngtab not found: %s", args.syngtab)
        return 2

    os.makedirs(args.out_dir, exist_ok=True)

    # Read syngtab into memory (header + rows)
    df = pd.read_csv(args.syngtab, sep="\t")
    
    out_path = os.path.join(args.out_dir, "busco_order.tsv")
    df['#marker'].to_csv(out_path, sep='\t', index=False, header=False)

    # For each node, write out file like the awk does
    for header in list(df.columns.values):
        if (header.startswith('n') and header.endswith("_seq")) is False:
            continue
        # print(header)
        out_path = os.path.join(args.out_dir, f"{header}.tsv")
        df[header].to_csv(out_path, sep='\t', index=False, header=False)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
