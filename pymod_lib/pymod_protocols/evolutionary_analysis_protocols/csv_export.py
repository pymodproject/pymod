import os
import csv
import pandas as pd
import warnings
#############################MODIFIED on  24/02/2025#################################
from PyQt5.QtWidgets import QFileDialog

from PyQt5.QtWidgets import QMessageBox




class CAMPOCSVExporter:
    def __init__(self, pymod, input_cluster_element, campo_window):
        """
        Initializes the CAMPOCSVExporter class with references to PyMod, 
        the input cluster element, and the CAMPO window.
        """
        self.pymod = pymod
        self.input_cluster_element = input_cluster_element
        self.campo_window = campo_window

    def save_csv(self, campo_list, selected_seq="All", selected_format="with alignment"):
        """
        Generates and saves the CAMPO scores as a CSV file, either for all sequences 
        or for a single selected sequence.
        """
        # Retrieve aligned sequences
        aligned_sequences = self.input_cluster_element.get_children()
        sequence_dict = {}  # Dictionary to store sequence data
        headers = ["Alignment Position"]  # Header for CSV file with "with alignment" format 

        # Populate the dictionary with sequence names and sequences
        for seq in aligned_sequences:
            seq_name = seq.my_header  # Sequence name (header)
            seq_sequence = seq.my_sequence  # Sequence itself
            sequence_dict[seq_name] = seq_sequence
            headers.append(seq_name)  # Add sequence name to header

        # Determine the maximum sequence length
        max_len = max(len(seq) for seq in sequence_dict.values())

        table = []  # List to store CSV rows

        ####################COMBINATION of sequence and format choices###############################

        #############################MODIFIED on  24/02/2025#################################
        #1. Check if the user selected 'No' (no export)
        if self.campo_window.save_file_cbx.get() == "No":
            # Do not generate CSV, just run the analysis
            print("CAMPO analysis completed without exporting results to CSV.")
            return  # Exit the function without generating CSV
        
        #2. If the user doesn't check "No"
        else:
            ###################################MODIFIED on 26/02/2025###########################################################
            #2.1 Check if the user selected 'All' sequences  and format with alignment
            if self.campo_window.save_file_cbx.get() == "All" and selected_format == "with alignment":

                #############################MODIFIED on  24/02/2025#################################
                # Ask the user for a directory and file name to save the CSV file
                save_dir = QFileDialog.getExistingDirectory(self.campo_window, "Select Directory")
                campo_txt_path = os.path.join(save_dir, f"alignment_table_campo_{selected_seq}_campo.csv")
        
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

            ###################################MODIFIED on 28/02/2025###########################################################
            ##2.2 The 'All' option supports  'without alignment' and save automatically a single csv file for each sequence of the alignments without the alignemnt information 
            elif self.campo_window.save_file_cbx.get() == "All" and selected_format == "without alignment":
                #Ask the user for the destination directory only once.
                save_dir = QFileDialog.getExistingDirectory(self.campo_window, "Select Directory")
                #Check that the user has selected a directory.
                if save_dir:  
                    #loop through each sequence in the sequence dict
                    for selected_seq in sequence_dict:
                        seq_sequence = sequence_dict[selected_seq]
                        table = []
                        
                        position = 1  # Start residue position from 1
                        
                        # Loop through the selected sequence
                        for i in range(len(seq_sequence)):
                            aa = seq_sequence[i]  # Residue at this position
                            score = campo_list[i]['campo-score'] if i < len(campo_list) else ""  # CAMPO score
                            
                            if aa != "-":
                                residue_position = f"AA{position}"  # Correct residue position
                                table.append([residue_position, f"{aa}", f"{score}"])
                                position += 1  # Increment only for valid residues
            
                        # 3 columns 
                        headers = [f"{selected_seq} Residue Position", "Residue", "CAMPO score"]
                        
                        #############################MODIFIED on  24/02/2025#################################
                        campo_txt_path = os.path.join(save_dir, f"{selected_seq}_campo_without_alignment_option_All.csv")
                        # Write the CSV file with the correct headers
                        with open(campo_txt_path, "w", newline="") as f:
                            writer = csv.writer(f)
                            writer.writerow(headers)
                            writer.writerows(table)
            

            ###################################MODIFIED on 26/02/2025###########################################################
            ##2.3 The 'selected_seq' option with 'with alignment' Format 
            elif self.campo_window.save_file_cbx.get() != "All" and selected_format == "with alignment":
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
                    
                    #############################MODIFIED on  24/02/2025#################################
                    # Ask the user for a directory and file name to save the CSV file
                    save_dir = QFileDialog.getExistingDirectory(self.campo_window, "Select Directory")
                    campo_txt_path = os.path.join(save_dir, f"{selected_seq}_campo_with_alignment.csv")
            
            ###################################MODIFIED on 26/02/2025###########################################################
            ##2.4 The 'selected_seq' option with 'without alignment' Format 
            #Code for saving csv file for a particular selected_seq without alignment information, only keeping the list of AA and their CAMPO score in the alignment 
            elif self.campo_window.save_file_cbx.get() != "All" and selected_format == "without alignment":

                if selected_seq in sequence_dict:
                    seq_sequence = sequence_dict[selected_seq]
                    table = []
                    
                    position = 1  # Start residue position from 1
                    
                    # Loop through the selected sequence
                    for i in range(len(seq_sequence)):
                        aa = seq_sequence[i]  # Residue at this position
                        score = campo_list[i]['campo-score'] if i < len(campo_list) else ""  # CAMPO score
                        
                        if aa != "-":
                            residue_position = f"AA{position}"  # Correct residue position
                            table.append([residue_position, f"{aa}", f"{score}"])
                            position += 1  # Increment only for valid residues
        
                    # 3 columns 
                    headers = [f"{selected_seq} Residue Position", "Residue", "CAMPO score"]
                    
                    #############################MODIFIED on  24/02/2025#################################
                    # Ask the user for a directory and file name to save the CSV file
                    save_dir = QFileDialog.getExistingDirectory(self.campo_window, "Select Directory")
                    campo_txt_path = os.path.join(save_dir, f"{selected_seq}_campo_without_alignment.csv")


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
