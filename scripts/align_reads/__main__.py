from .read_aligner import ReadAligner

def main():
    # Create an instance of the ReadAligner class
    read_aligner = ReadAligner()

    # Read input data
    read_aligner.read_inputs()

    # Map clusters
    read_aligner.map_clusters()


if __name__ == '__main__':
    main()
