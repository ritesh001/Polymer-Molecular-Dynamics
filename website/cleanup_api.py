import glob
import shutil


def cleanup(cleanup_route: str):
    for file in glob.glob(f'{cleanup_route}/*.md'):
        data = []
        # read input file
        with open(file, 'rt') as f:
            for line in f:
                # replace all occurrences of the required string
                line = line.replace('####', '###')
                line = line.replace('&quot;', '"')
                line = line.replace('&lt;', '<')
                data.append(line)

        # open the input file in write mode
        with open(file, 'wt') as f:
            f.writelines(data)


if __name__ == '__main__':
    # clean up the main modules
    cleanup('./api/core')
    cleanup('./api/postprocessing')

    # don't show the util module in the doc
    try:
        shutil.rmtree('./api/util')
    except Exception:
        pass
