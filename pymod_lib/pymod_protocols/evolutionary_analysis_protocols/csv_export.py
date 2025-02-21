import os
import csv
import pandas as pd
import warnings

class CAMPOCSVExporter:
    def __init__(self, pymod, input_cluster_element, campo_window):
        """
        Initializes the CAMPOCSVExporter class with references to PyMod, 
        the input cluster element, and the CAMPO window.
        """
        self.pymod = pymod
        self.input_cluster_element = input_cluster_element
        self.campo_window = campo_window

    def save_csv(self, campo_list, selected_seq="All"):
        """
        Generates and saves the CAMPO scores as a CSV file, either for all sequences 
        or for a single selected sequence.
        """
        # Retrieve aligned sequences
        aligned_sequences = self.input_cluster_element.get_children()
        sequence_dict = {}  # Dictionary to store sequence data
        headers = ["Alignment Position"]  # Header for CSV file

        # Populate the dictionary with sequence names and sequences
        for seq in aligned_sequences:
            seq_name = seq.my_header  # Sequence name (header)
            seq_sequence = seq.my_sequence  # Sequence itself
            sequence_dict[seq_name] = seq_sequence
            headers.append(seq_name)  # Add sequence name to header

        # Determine the maximum sequence length
        max_len = max(len(seq) for seq in sequence_dict.values())

        table = []  # List to store CSV rows

        # Check if the user selected 'All' sequences or a specific one
        if self.campo_window.save_file_cbx.get() == "All":
            # File path for saving all sequences
            campo_txt_path = os.path.join(self.pymod.alignments_dirpath, f"alignment_table_campo_{selected_seq}.csv")
            
            # Loop through each alignment position
            for i in range(max_len):
                row = [f"AA{i+1}"]  # Alignment position (AA1, AA2, ...)
                
                for seq_id in sequence_dict:
                    seq = sequence_dict[seq_id]
                    if i < len(seq):
                        aa = seq[i]  # Get amino acid at position i
                        score = campo_list[i]['campo-score'] if i < len(campo_list) else ""
                        row.append(f"{aa},{score}" if score else f"{aa},")  # Append score if available
                    else:
                        row.append("")  # Leave empty if sequence is shorter

                table.append(row)  # Add row to table
        
        else:
            # If a specific sequence is selected
            if selected_seq in sequence_dict:
                seq_sequence = sequence_dict[selected_seq]
                table = []
                
                # Loop through the selected sequence
                for i in range(len(seq_sequence)):
                    alignment_position = f"AA{i+1}"  # Alignment position
                    aa = seq_sequence[i]  # Residue at this position
                    score = campo_list[i]['campo-score'] if i < len(campo_list) else ""

                    if aa == "-":
                        table.append([alignment_position, f"{aa}"])  # No score for gaps
                    else:
                        table.append([alignment_position, f"{aa},{score}" if score else f"{aa},"])

                # Update header to include only the selected sequence
                headers = ["Alignment Position", selected_seq]
                
                # File path for saving the selected sequence
                campo_txt_path = os.path.join(self.pymod.alignments_dirpath, f"{selected_seq}_campo.csv")

        # Write the CSV file with the correct headers
        with open(campo_txt_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(headers)
            writer.writerows(table)

        # Format the CSV file correctly
        self.right_format_csv(campo_txt_path)

    def right_format_csv(self, table_campo):
        """
        Formats the CSV file to correct numbering and ensure consistency.
        """
        #Overwrite the CSV file with the corrected format
        try:
            df = pd.read_csv(table_campo)  # Read the CSV file

            # Process each column except the first (alignment position)
            for col in df.columns[1:]:
                counter = 1  # Counter for numbering amino acids
                for i in range(len(df)):
                    if isinstance(df[col][i], str) and "," in df[col][i]:
                        aa, score = df[col][i].split(",")  # Split amino acid and score
                        try:
                            score = float(score)  # Convert score to float
                        except ValueError:
                            score = 0.0  # Fallback if conversion fails

                        if aa == "-":
                            df.at[i, col] = "-"  # Keep gaps unchanged
                        else:
                            df.at[i, col] = f"{aa}-{counter},{score}"  # Format as 'AA-Position,Score'
                            counter += 1  # Increment position counter

            df.to_csv(table_campo, index=False)  # Save formatted CSV file
        
        #warnings if the file doens't exist or no permissions to write it or other error
        except FileNotFoundError:
            warnings.warn(f"File not found: {table_campo}.", RuntimeWarning)
        except PermissionError:
            warnings.warn(f"Permission denied: Cannot write to {table_campo}.", RuntimeWarning)
        except Exception as e:
            warnings.warn(f"An error occurred: {e}", RuntimeWarning)
