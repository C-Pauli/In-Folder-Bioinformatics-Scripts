import os
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

def collapse_exons(features):
    """Collapse exons into single genes, separating genes if intron length > 5000 bp."""
    collapsed_features = []
    current_gene = None
    current_gene_start = None
    current_gene_end = None
    
    for feature in features:
        if feature.type == "CDS":
            if current_gene is None:
                # Start a new gene
                current_gene = feature
                current_gene_start = feature.location.start
                current_gene_end = feature.location.end
            else:
                # If the next CDS is part of the same gene (within 5000 bp), extend the current gene
                intron_length = feature.location.start - current_gene_end
                if intron_length <= 5000:
                    current_gene_end = feature.location.end
                else:
                    # Save the current gene and start a new one
                    current_gene.location = FeatureLocation(current_gene_start, current_gene_end)
                    collapsed_features.append(current_gene)
                    current_gene = feature
                    current_gene_start = feature.location.start
                    current_gene_end = feature.location.end
    
    # Don't forget the last gene
    if current_gene is not None:
        current_gene.location = FeatureLocation(current_gene_start, current_gene_end)
        collapsed_features.append(current_gene)

    return collapsed_features

def process_gbk_files():
    """Process all .gbk files in the current directory and output the modified files."""
    # Get the current working directory
    current_dir = os.getcwd()

    # Get all .gbk files in the current directory
    gbk_files = [f for f in os.listdir(current_dir) if f.endswith(".gbk")]

    for gbk_file in gbk_files:
        print(f"Processing {gbk_file}...")
        
        input_file = os.path.join(current_dir, gbk_file)
        output_file = os.path.join(current_dir, gbk_file.replace(".gbk", "_processed.gbk"))
        
        # Parse the input GenBank file
        with open(input_file, "r") as input_handle:
            records = SeqIO.parse(input_handle, "genbank")
            
            processed_records = []
            for record in records:
                # Process each record (each sequence in the file)
                collapsed_features = collapse_exons(record.features)
                # Create a new record with collapsed features
                new_record = record
                new_record.features = collapsed_features
                processed_records.append(new_record)

            # Write the processed records to a new GenBank file
            with open(output_file, "w") as output_handle:
                SeqIO.write(processed_records, output_handle, "genbank")
        
        print(f"Finished processing {gbk_file}, output saved to {output_file}")

# Process the .gbk files in the current directory
process_gbk_files()

