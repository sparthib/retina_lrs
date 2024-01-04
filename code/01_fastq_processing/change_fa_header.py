with open('/dcs04/hicks/data/sparthib/references/sample.fa', 'r') as file:
    for line in file:
        if line.startswith('>'):
            header_fields = line.strip().split('|')
            first_field = header_fields[0][1:]  # Extracting the first field after '>'
            updated_lines.append(f'>{first_field}\n')  # Creating the updated header line
        else:
            updated_lines.append(line)

# Write the updated content back to the file
with open(file_path, 'w') as file:
    file.writelines(updated_lines)
