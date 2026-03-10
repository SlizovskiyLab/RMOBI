import csv

csv_file_path = "../data/MGE_total_classification.xlsx - Sheet.csv"
text_file_path = '../scripts/mge_id.txt'
output_file_name = '../scripts/mge_name.txt'

'''
for plasmids/prophages (where final_classification= plasmid), using 'socus to name the gene
for both plasmids / prophages, if the final_classificaiton= likelyIS/TE, then the gene named based on sig_seq column
For ICEBerg or anything in the final_classification= ICE, using the genebank_id
For virus, classify  by gene:vir 
'''


try:
    csv_data = {}
    with open(csv_file_path, mode='r', newline='', encoding='utf-8') as infile:
        reader = csv.reader(infile)
        header = next(reader)
        try:
            id_col_index = header.index('IDs')
            final_class_col_index = header.index('final_classification')
            sig_seq_col_index = header.index('sig_seq')
            socus_col_index = header.index('socus')
            genebank_id_index = header.index('genebank_id')
        except ValueError as e:
            print(f"Error: A required column is missing from the CSV header: {e}")
            exit()

        for row in reader:
            if len(row) > max(id_col_index, final_class_col_index, sig_seq_col_index, genebank_id_index, socus_col_index):
                key = row[id_col_index]
                final_classification = row[final_class_col_index]
                sig_seq = row[sig_seq_col_index]
                genebank_id = row[genebank_id_index]
                socus = row[socus_col_index]

                value_to_store = None

                if final_classification == "likely IS/TE":
                    value_to_store = sig_seq
                elif final_classification == "ICE":
                    value_to_store = genebank_id
                elif socus and socus.strip():
                    value_to_store = socus

                if value_to_store and value_to_store.strip():
                    csv_data[key] = value_to_store

    with open(text_file_path, mode='r') as infile, open(output_file_name, mode='w') as outfile:
        for line in infile:
            cleaned_line = line.strip().rstrip(',').lstrip('{').rstrip('}')
            parts = cleaned_line.split(',', 1)

            if len(parts) == 2:
                number = parts[0].strip()
                key_from_text = parts[1].strip().strip('"')

                if key_from_text in csv_data:
                    new_value = csv_data[key_from_text]
                    output_line = f'{{{number}, "{new_value}"}},'
                    outfile.write(output_line + '\n')
                else:
                    outfile.write(line)
            else:
                outfile.write(line)

    print(f"Processing complete. Data saved to {output_file_name}")

except FileNotFoundError as e:
    print(f"Error: One of the files was not found. Please check the paths: {e}")
except Exception as e:
    print(f"An error occurred: {e}")