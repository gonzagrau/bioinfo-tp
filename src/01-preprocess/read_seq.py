from Bio import SeqIO
import pandas as pd
import os


def extract_cds_features(record):
    """
    Extract CDS (coding sequence) features from a GenBank record.

    Args:
        record: A SeqIO record object from a GenBank file

    Returns:
        A list of dictionaries containing CDS information
    """
    cds_features = []

    for feature in record.features:
        if feature.type == "CDS":
            # Handle split locations (introns/exons)
            if hasattr(feature.location, "parts"):
                locations = []
                for part in feature.location.parts:
                    locations.append(f"{part.start}..{part.end}")
                location_str = f"join({','.join(locations)})"
            else:
                location_str = f"{feature.location.start}..{feature.location.end}"

            # Extract metadata
            cds_seq = feature.location.extract(record.seq)
            qualifiers = feature.qualifiers
            gene_name = qualifiers.get("gene", ["Unknown"])[0]
            product = qualifiers.get("product", ["Unknown"])[0]
            protein_id = qualifiers.get("protein_id", ["Unknown"])[0]

            # Get translation, ideally from qualifiers, if available
            if "translation" in qualifiers:
                translation = qualifiers["translation"][0]
            else:
                translation = cds_seq.translate()

            cds_features.append(
                {
                    "gene": gene_name,
                    "location": location_str,
                    "product": product,
                    "protein_id": protein_id,
                    "cds_sequence": str(cds_seq),
                    "translation": str(translation),
                }
            )

    return cds_features


def save_cds_data(cds_features, record_name, output_dir="../../outputs/fastas"):
    """
    Processes and saves coding sequence (CDS) features and their associated data for a specified record into files.
    It creates two FASTA files containing both CDS sequences and their translations, and also a CSV file with CDS
    feature details. The files are saved in the specified output directory.

    :param cds_features: List of dictionaries, where each dictionary represents a CDS feature. Each dictionary should
                         contain keys such as 'gene', 'protein_id', 'product', 'cds_sequence', and 'translation'.
    :param record_name: The name of the record used as a basis for naming the output files.
    :param output_dir: The directory where the output files will be saved. Defaults to 'output'.
    :type output_dir: str
    :return: None
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    df = pd.DataFrame(cds_features)
    df.to_csv(f"{output_dir}/{record_name}_cds_features.csv", index=False)

    # now: one file per CDS *and* one per protein
    for i, feature in enumerate(cds_features, start=1):
        # build a safe filename: include record, index or protein_id
        base = f"{record_name}_{feature['protein_id']}"
        
        # write the CDS sequence
        cds_path = os.path.join(output_dir, f"{base}_cds.fasta")
        with open(cds_path, "w") as f_cds:
            header = f">{feature['gene']}|{feature['protein_id']}|{feature['product']}"
            f_cds.write(header + "\n")
            f_cds.write(feature['cds_sequence'] + "\n")

        # write the translated protein
        prot_path = os.path.join(output_dir, f"{base}_protein.fasta")
        with open(prot_path, "w") as f_prot:
            f_prot.write(header + "\n")
            f_prot.write(feature['translation'] + "\n")


def make_joinded_files(output_dir="../../outputs/fastas"):
    """
    Joins all CDS and protein FASTA files into single files for each type.
    This function assumes that the output directory contains the individual CDS and protein files.

    :param output_dir: The directory where the output files are saved. Defaults to 'output'.
    :type output_dir: str
    :return: None
    """
    cds_files = [f for f in os.listdir(output_dir) if f.endswith("_cds.fasta")]
    prot_files = [f for f in os.listdir(output_dir) if f.endswith("_protein.fasta")]
    joined_dir = os.path.join(output_dir, "joined")
    if not os.path.exists(joined_dir):
        os.makedirs(joined_dir)

    with open(os.path.join(joined_dir, "all_cds.fasta"), "w") as f_cds_all:
        for cds_file in cds_files:
            with open(os.path.join(output_dir, cds_file)) as f_cds:
                f_cds_all.write(f_cds.read())

    with open(os.path.join(joined_dir, "all_proteins.fasta"), "w") as f_prot_all:
        for prot_file in prot_files:
            with open(os.path.join(output_dir, prot_file)) as f_prot:
                f_prot_all.write(f_prot.read())


def main():
    for file in os.listdir("../../data"):
        # Skip irrelevant files
        if not file.endswith(".gb") or file.endswith(".gbk"):
            continue

        record_name = file.split(".")[0]
        print(f"Reading GenBank file: {record_name}")
        with open(os.path.join("../../data", file)) as handle:
            record = SeqIO.read(handle, format="genbank")

        cds_features = extract_cds_features(record)
        print(f"Found {len(cds_features)} CDS features in the GenBank file:\n")

        for i, feature in enumerate(cds_features):
            print(f"CDS {i+1}:")
            print(f"  Gene: {feature['gene']}")
            print(f"  Product: {feature['product']}")
            print(f"  Location: {feature['location']}")
            print(f"  Protein ID: {feature['protein_id']}")
            print(f"  CDS Sequence: {feature['cds_sequence'][:50]}...")
            print(f"  Translation: {feature['translation'][:50]}...")
            print()

        # Save data to files
        save_cds_data(cds_features, record_name)
        print(f"CDS data saved to the 'output' directory.")

    # Join all CDS and protein files into single files
    make_joinded_files(output_dir="../../outputs/fastas")
    print("Joined CDS and protein files created in the 'joined' directory.")




if __name__ == "__main__":
    main()
