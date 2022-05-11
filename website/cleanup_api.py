import glob


def cleanup(cleanup_route: str):
    for file in glob.glob(f'{cleanup_route}/*.md'):
        #read input file
        fin = open(file, 'rt')
        #read file contents to string
        data = fin.read()
        #replace all occurrences of the required string
        data = data.replace('####', '###')
        data = data.replace('&quot;', '"')
        data = data.replace('&lt;', '<')
        #close the input file
        fin.close()
        #open the input file in write mode
        fin = open(file, 'wt')
        #overrite the input file with the resulting data
        fin.write(data)
        #close the file
        fin.close()


if __name__ == '__main__':
    cleanup('./api/core')
    cleanup('./api/postprocessing')
