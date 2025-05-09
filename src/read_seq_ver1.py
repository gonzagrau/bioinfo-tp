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


def save_cds_data(cds_features, record_name, output_dir="../output/fastas"):
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

    # Save sequences to FASTA files
    with open(f"{output_dir}/{record_name}_cds_sequences.fasta", "w") as f:
        for i, feature in enumerate(cds_features):
            f.write(
                f">{feature['gene']}|{feature['protein_id']}|{feature['product']}\n"
            )
            f.write(f"{feature['cds_sequence']}\n")

    # Save translations to FASTA files
    with open(f"{output_dir}/{record_name}_protein_sequences.fasta", "w") as f:
        for i, feature in enumerate(cds_features):
            f.write(
                f">{feature['gene']}|{feature['protein_id']}|{feature['product']}\n"
            )
            f.write(f"{feature['translation']}\n")


def main():
    for file in os.listdir("../data"):
        # Skip irrelevant files
        if not file.endswith(".gb") or file.endswith(".gbk"):
            continue

        record_name = file.split(".")[0]
        print(f"Reading GenBank file: {record_name}")
        with open(os.path.join("../data", file)) as handle:
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


if __name__ == "__main__":
    main()
