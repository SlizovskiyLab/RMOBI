import csv

csv_file_path = "../../data/MGE_total_classification.xlsx - Sheet.csv"
text_file_path = '../../scripts/mge_id.txt'
output_file_name = '../../scripts/mge_group.txt'

try:
    csv_data = {}
    with open(csv_file_path, mode='r', newline='') as infile:
        reader = csv.reader(infile)
        next(reader)
        for row in reader:
            if len(row) >= 2:
                csv_data[row[0]] = row[1]

    with open(text_file_path, mode='r') as infile, open(output_file_name, mode='w') as outfile:
        for line in infile:
            cleaned_line = line.strip().rstrip(',').lstrip('{').rstrip('}')
            parts = cleaned_line.split(',', 1)
            
            if len(parts) == 2:
                number = parts[0].strip()
                key = parts[1].strip().strip('"')

                if key in csv_data:
                    new_value = csv_data[key]
                    output_line = f'{{{number}, "{new_value}"}},'
                    outfile.write(output_line + '\n')
                else:
                    outfile.write(line)

    print(f"Data saved to {output_file_name}")

except FileNotFoundError:
    print(f"Error: One of the files was not found.")
except Exception as e:
    print(f"An error occurred: {e}")