# Project Structure

'''
OPropagator/
├── docs/                       # Documentation
│   └── overview.md             # Explanation of the model and its features
├── orbit/                      # Main package for core functionality
│   ├── __init__.py             # Initialize package
│   ├── forces.py               # Gravitational and other forces (e.g., atmospheric drag)
│   ├── initial_conditions.py   # Functions to define or read initial conditions
│   ├── models.py               # Orbital models (e.g., two-body, n-body)
│   ├── propagator.py           # Core propagation functions (e.g., RK4, Euler)
│   └── utils.py                # Helper functions (e.g., unit conversions)
├── scripts/                    # Scripts to run specific tasks (e.g., examples)
│   ├── n_body_example.py       # Example of n-body simulation
│   └── single_orbit_example.py # Example of single orbit propagation
├── tests/                      # Test suite
│   ├── __init__.py
│   ├── test_forces.py          # Tests for forces functions
│   ├── test_models.py          # Tests for orbital models
│   └── test_propagator.py      # Tests for propagator functions
├── .gitignore                  # Ignore unnecessary files in version control
├── LICENSE.md                  # License
├── main.py                     # Main entry point of the program
├── README.md                   # Project overview and setup instructions
└── requirements.txt            # List of dependencies
'''