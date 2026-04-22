from .read_aligner import ReadAligner

def main():
    # Create an instance of the ReadAligner class
    read_aligner = ReadAligner()

    # Read input data
    read_aligner.read_inputs()

    # Map clusters
    # read_aligner.map_clusters()
    
    # MOUMI: added option to map a single cluster for parallelization
    if read_aligner.inputs.get("single_cluster", False):
        read_aligner.map_single_cluster()
    else:
        read_aligner.map_clusters()


if __name__ == '__main__':
    main()
