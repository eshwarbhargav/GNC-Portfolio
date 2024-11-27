from src.orbit.propagator import FakeFixedStepHandler

# Set up the Propagator
propagator_test = [propagator_eh, propagator_num]

x_positions = []
y_positions = []
z_positions = []

for propagator in propagator_test:
    # Set up the Handler
    handler = FakeFixedStepHandler()
    propagator.getMultiplexer().add(initStep, handler)

    # Propagate for 24 hours
    finalState = propagator.propagate(initialDate, initialDate.shiftedBy(7*24*60.0*60.0))
    # X, Y, Z positions in meters
    x_positions.append(handler.x_positions)    # List of X positions from propagation
    y_positions.append(handler.y_positions)    # List of Y positions from propagation
    z_positions.append(handler.z_positions)    # List of Z positions from propagation

tpva_start = initialState.getPVCoordinates().getPosition()    # Initial position
tpva_end = finalState.getPVCoordinates().getPosition()    # Final position
