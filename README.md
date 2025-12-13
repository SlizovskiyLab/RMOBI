# RMOBI: a Network Framework for Tracking Resistome-Mobilome Colocalization Across Fecal Microbiota Transplantation Cohorts

RMOBI is a framework built in C++ for constructing and analyzing large scale heterogeneous temporal networks of antimicrobial resistance gene (ARG) and mobile genetic element (MGE) colocalizations across fecal microbiota transplantation (FMT) phases. The tool models colocalization and temporal relationships from patient metagenomic data, producing structured outputs for interactive visualization using D3.js.

## Project Overview

The goal of RMOBI is to model and visualize the dynamic interactions between ARGs and MGEs before and after FMT treatment. The C++ implementation generates structured JSON files containing node and edge information, representing both colocalization edges and temporal edges that capture ARG–MGE evolution across donor, pre-FMT, and post-FMT phases.

The accompanying web based visualization module, built with D3.js, dynamically renders these networks, enabling users to explore treatment related trends and filter the network interactively based on clinical and biological criteria.


## Features

1. Efficient C++ implementation for constructing heterogeneous temporal graphs.
   
2.. Supports multiple categorization modes, including disease type and temporal dynamics.

3. Produces structured CSV files that summarize Disease specific, MGE group and temporal dynamics of colocalizations before and after FMT.

4. Interactive D3.js visualization with real time filtering and multi phase exploration.

5. Filter by disease type, patient phase, MGE group, or drill down to individual ARG/MGE nodes.

6. Export network visualizations as SVG files for publication or downstream analysis.

7. Fully documented workflow and configuration examples for reproducibility.

---

### The web-based system will be available on the following link

## Run the Program  

### Prerequisites

C++17 or later

Make (for build automation)

Standard C++ libraries (<iostream>, <fstream>, <map>, <vector>, etc.)

JSON library (e.g., nlohmann/json)


# Clone the repository
git clone https://github.com/<your-repo>/AMRGraphNet.git
cd AMRGraphNet

# Build the C++ project
make



### Linux/macOS
Step 1: Create and enter build directory
mkdir -p build
cd build

Step 2: Generate Makefiles with CMake
```
cmake ..
```

Step 3: Compile the project
```
cmake --build .
```

Step 4: Run the executable
```
./CoNet
```

### Windows
Step 1: Create and enter build directory
mkdir -p build
cd build

Step 2: Generate Makefiles with CMake
```
cmake .. -G "MinGW Makefiles"
```

Step 3: Compile the project
```
mingw32-make
```

Step 4: Run the executable
```
./CoNet.exe
```

### Manually Run
g++ -std=c++17 -Wall -Iinclude -Ithird_party src/*.cpp -o CoNet.exe

./CoNet.exe


cd viz
python3 -m http.server 8080

Check 
http://localhost:8080/index.html




