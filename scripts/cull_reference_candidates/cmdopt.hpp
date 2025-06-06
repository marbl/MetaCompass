#include <vector>

using namespace std;

class Cmdopt {
    public:
        string reads, output_folder, reference_folders;
        int min_shared_kmer;
        Cmdopt();
        static void refCullingHelp(int status);
        void parseCmdOptions(int argc, char* argv[], Cmdopt* cmdopt);
};