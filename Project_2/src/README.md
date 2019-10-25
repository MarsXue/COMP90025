## COMP90025 Project 2
The focus of this project is to simulate the two-dimensional N-body problem.

### Team Member
- Hanbin Li - $hanbinl1$
- Wenqing Xue - $wenqingx$

### Implementation
- Basic OpenMP + MPI approach: $N^2$ complexity
  - N2_code.c
- Barnes-Hut OpenMP + MPI approach: $N\log N$ complexity
  - NlogN_code.c
  - tree.c

### Folder Structure
```js
/
├── figure.py           // Python code for figures
├── generator.c         // C code for generating the test data
├── N2/                 // Folder for N^2 approach
├── NlogN/              // Folder for Barnes-Hut approach
├── test_1200.txt       // Test data with 1200 particles
├── test_2400.txt       // Test data with 2400 particles
├── test_3600.txt       // Test data with 3600 particles
├── test_4800.txt       // Test data with 4800 particles
└── test_6000.txt       // Test data with 6000 particles
```