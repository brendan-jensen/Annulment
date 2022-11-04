import re


class HeaderReplacer:
    def __init__(self, filename):
        self.filename = filename

    def replace_header(self, indices):
        with open(self.filename, 'r') as file:
            file_data = file.read()

        header = re.findall("CSD_code(.*)", file_data)[0]

        replace_string = "CSD_code" + header
        # new_string = str(index)
        new_string = ','.join([str(i) for i in indices])

        file_data = file_data.replace(replace_string, new_string)

        with open(self.filename, 'w') as file:
            file.write(file_data)
