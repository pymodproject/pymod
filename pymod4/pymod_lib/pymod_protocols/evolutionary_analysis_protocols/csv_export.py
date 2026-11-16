#Created by SARA GIORDANI

import os
import csv
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
        Generates and saves the CAMPO scores as CSV file(s), then shows a
        confirmation pop-up listing every file that was written.
        Returns the list of saved file paths (empty if nothing was saved).
        """
        aligned_sequences = self.input_cluster_element.get_children()
        sequence_dict = {}
        headers = ["Alignment Position"]

        for seq in aligned_sequences:
            sequence_dict[seq.my_header] = seq.my_sequence
            headers.append(seq.my_header)

        max_len = max(len(seq) for seq in sequence_dict.values())

        table = []
        saved_files = []  # every CSV actually written ends up here

        # 1. User chose "No": run analysis, export nothing.
        if self.campo_window.save_file_cbx.get() == "No":
            print("CAMPO analysis completed without exporting results to CSV.")
            return saved_files

        # 2.1 All + with alignment (single combined file).
        if self.campo_window.save_file_cbx.get() == "All" and selected_format == "with alignment":
            save_dir = QFileDialog.getExistingDirectory(self.campo_window, "Select Directory")
            if not save_dir:                      # user cancelled the dialog
                return saved_files
            campo_txt_path = os.path.join(save_dir, f"alignment_table_campo_{selected_seq}_campo.csv")

            for i in range(max_len):
                row = [f"AA{i+1}"]
                for seq_id in sequence_dict:
                    seq = sequence_dict[seq_id]
                    if i < len(seq):
                        aa = seq[i]
                        score = campo_list[i]['campo-score'] if i < len(campo_list) else ""
                        row.append(f"{aa},{score}" if score else f"{aa},")
                    else:
                        row.append("")
                table.append(row)

            self._write_csv(campo_txt_path, headers, table)
            self.right_format_csv(campo_txt_path)
            saved_files.append(campo_txt_path)

        # 2.2 All + without alignment (one file per sequence).
        elif self.campo_window.save_file_cbx.get() == "All" and selected_format == "without alignment":
            save_dir = QFileDialog.getExistingDirectory(self.campo_window, "Select Directory")
            if not save_dir:
                return saved_files
            for seq_name in sequence_dict:
                seq_sequence = sequence_dict[seq_name]
                table = []
                position = 1
                for i in range(len(seq_sequence)):
                    aa = seq_sequence[i]
                    score = campo_list[i]['campo-score'] if i < len(campo_list) else ""
                    if aa != "-":
                        table.append([f"AA{position}", f"{aa}", f"{score}"])
                        position += 1
                file_headers = [f"{seq_name} Residue Position", "Residue", "CAMPO score"]
                campo_txt_path = os.path.join(save_dir, f"{seq_name}_campo_without_alignment_option_All.csv")
                self._write_csv(campo_txt_path, file_headers, table)
                saved_files.append(campo_txt_path)

        # 2.3 selected_seq + with alignment.
        elif self.campo_window.save_file_cbx.get() != "All" and selected_format == "with alignment":
            if selected_seq in sequence_dict:
                seq_sequence = sequence_dict[selected_seq]
                table = []
                for i in range(len(seq_sequence)):
                    alignment_position = f"AA{i+1}"
                    aa = seq_sequence[i]
                    score = campo_list[i]['campo-score'] if i < len(campo_list) else ""
                    if aa == "-":
                        table.append([alignment_position, f"{aa}"])
                    else:
                        table.append([alignment_position, f"{aa},{score}" if score else f"{aa},"])
                file_headers = ["Alignment Position", selected_seq]
                save_dir = QFileDialog.getExistingDirectory(self.campo_window, "Select Directory")
                if not save_dir:
                    return saved_files
                campo_txt_path = os.path.join(save_dir, f"{selected_seq}_campo_with_alignment.csv")
                self._write_csv(campo_txt_path, file_headers, table)
                self.right_format_csv(campo_txt_path)
                saved_files.append(campo_txt_path)

        # 2.4 selected_seq + without alignment.
        elif self.campo_window.save_file_cbx.get() != "All" and selected_format == "without alignment":
            if selected_seq in sequence_dict:
                seq_sequence = sequence_dict[selected_seq]
                table = []
                position = 1
                for i in range(len(seq_sequence)):
                    aa = seq_sequence[i]
                    score = campo_list[i]['campo-score'] if i < len(campo_list) else ""
                    if aa != "-":
                        table.append([f"AA{position}", f"{aa}", f"{score}"])
                        position += 1
                file_headers = [f"{selected_seq} Residue Position", "Residue", "CAMPO score"]
                save_dir = QFileDialog.getExistingDirectory(self.campo_window, "Select Directory")
                if not save_dir:
                    return saved_files
                campo_txt_path = os.path.join(save_dir, f"{selected_seq}_campo_without_alignment.csv")
                self._write_csv(campo_txt_path, file_headers, table)
                saved_files.append(campo_txt_path)

        # Single confirmation pop-up listing everything that was written.
        self._show_save_confirmation(saved_files)
        return saved_files

    def _write_csv(self, path, headers, rows):
        """Writes a header row plus data rows to a CSV file."""
        with open(path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(headers)
            writer.writerows(rows)

    def _show_save_confirmation(self, saved_files):
        """Shows a pop-up confirming which CSV file(s) were saved."""
        if not saved_files:
            return  # nothing saved (e.g. user cancelled) -> no false confirmation
        if len(saved_files) == 1:
            message = "The CSV file was saved successfully:\n\n%s" % saved_files[0]
        else:
            message = ("The following %d CSV files were saved successfully:\n\n%s"
                       % (len(saved_files), "\n".join(saved_files)))
        QMessageBox.information(self.campo_window, "CSV Export Complete", message)

    def right_format_csv(self, table_campo):
        """
        Formats the CSV file to correct numbering and ensure consistency.
        """
        #Overwrite the CSV file with the corrected format
        try:
            # Read the CSV file (first row is the header).
            with open(table_campo, "r", newline="") as f:
                rows = list(csv.reader(f))

            if not rows:
                return  # Empty file, nothing to format.

            header = rows[0]
            data_rows = rows[1:]
            num_cols = len(header)

            # Process each column except the first (alignment position).
            for col in range(1, num_cols):
                counter = 1  # Counter for numbering amino acids
                for row in data_rows:
                    if col >= len(row):
                        continue  # Ragged row: this sequence is shorter.
                    cell = row[col]
                    if isinstance(cell, str) and "," in cell:
                        aa, score = cell.split(",")  # Split amino acid and score
                        try:
                            score = float(score)  # Convert score to float
                        except ValueError:
                            score = 0.0  # Fallback if conversion fails

                        if aa == "-":
                            row[col] = "-"  # Keep gaps unchanged
                        else:
                            row[col] = f"{aa}-{counter},{score}"  # Format as 'AA-Position,Score'
                            counter += 1  # Increment position counter

            # Save formatted CSV file.
            with open(table_campo, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(header)
                writer.writerows(data_rows)
        
        #warnings if the file doens't exist or no permissions to write it or other error
        except FileNotFoundError:
            warnings.warn(f"File not found: {table_campo}.", RuntimeWarning)
        except PermissionError:
            warnings.warn(f"Permission denied: Cannot write to {table_campo}.", RuntimeWarning)
        except Exception as e:
            warnings.warn(f"An error occurred: {e}", RuntimeWarning)
