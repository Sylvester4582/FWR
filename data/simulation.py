import pterasoftware as ps
import numpy as np
import matplotlib.pyplot as plt
from Selfpunctions import extract_forces, extract_second_cycle
import os 
import csv


def Mk4SaveDataToCSV(FlappingPeriod, Va, AoA, Lift, PitchingMoment, InducedDrag, file_name):
    # Ensure that Lift, PitchingMoment, and InducedDrag are numpy arrays
    Lift = np.asarray(Lift)
    PitchingMoment = np.asarray(PitchingMoment)
    InducedDrag = np.asarray(InducedDrag)

    # Total number of points for one flapping period
    num_points = Lift.size
    
    # Ensure that all arrays have the same length
    if not all(arr.size == num_points for arr in [PitchingMoment, InducedDrag]):
        raise ValueError(f"Lift, PitchingMoment, and InducedDrag must all have the same length of {num_points} points.")
    
    # Create a normalized time array (range between 0 and 1)
    normalized_time = np.linspace(0, 1, num_points, endpoint=False)
    
    # Write the CSV file in the current directory
    with open(file_name, mode='w', newline='') as file:
        writer = csv.writer(file)
        
        # Write the header (column names)
        writer.writerow(['Flapping Frequency', 'AirSpeed', 'Angle of Attack', 'Normalised Time', 'Lift', 'Induced Drag', 'Pitching Moment'])
        
        # Write the data for each time step
        for i in range(num_points):
            writer.writerow([FlappingPeriod, Va, AoA, normalized_time[i], Lift[i], InducedDrag[i], PitchingMoment[i]])

    print(f"Data has been successfully saved to {file_name}")

    
# Function call for the unsteady aerodynamic model with the NACA 8304 airfoil
def simulation(va= 5, aoa = 5, fp = 0.2):
    example_airplane=ps.geometry.Airplane(
        name="naca8304",
        x_ref=0.11,
        y_ref=0.0,
        z_ref=0.0,
        s_ref=None,
        b_ref=None,
        c_ref=None,
        wings=[
            ps.geometry.Wing(
                name="Main Wing",
                x_le=0.0,
                y_le=0.0,
                z_le=0.0,
                symmetric=True,
                num_chordwise_panels=6,
                chordwise_spacing="cosine",
                wing_cross_sections=[
                    ps.geometry.WingCrossSection(
                        x_le=0.0,
                        y_le=0.0,
                        z_le=0.0,
                        twist=0.0,
                        control_surface_type="symmetric",
                        control_surface_hinge_point=0.0,
                        control_surface_deflection=0.0,
                        num_spanwise_panels=8,
                        spanwise_spacing="cosine",
                        chord=0.3,
                        airfoil=ps.geometry.Airfoil(
                            name="naca8304",
                            coordinates=None,
                            repanel=True,
                            n_points_per_side=400,
                        ),
                    ),
                    ps.geometry.WingCrossSection(
                        x_le=0.0,
                        y_le=0.7,
                        z_le=0.0,
                        chord=0.2,
                        twist=0.0,
                        airfoil=ps.geometry.Airfoil(
                            name="naca8304",
                        ),
                    ),
                ],
            ),
            ps.geometry.Wing(
                name="V-Tail",
                x_le=0.45,
                y_le=0.0,
                z_le=0.0,
                num_chordwise_panels=6,
                chordwise_spacing="cosine",
                symmetric=True,
                wing_cross_sections=[
                    ps.geometry.WingCrossSection(
                        chord=0.2,
                        control_surface_type="symmetric",
                        control_surface_hinge_point=0.1,
                        control_surface_deflection=0.0,
                        airfoil=ps.geometry.Airfoil(
                            name="naca0004",
                        ),
                        twist=0.0,
                    ),
                    ps.geometry.WingCrossSection(
                        x_le=0.19,
                        y_le=0.2,
                        z_le=0.003,
                        chord=0.01,
                        control_surface_type="symmetric",
                        control_surface_hinge_point=0.1,
                        control_surface_deflection=0.0,
                        twist=0.0,
                        airfoil=ps.geometry.Airfoil(
                            name="naca0004",
                        ),
                    ),
                ],
            ),
            # Define the next wing.
            ps.geometry.Wing(
                name="V-Tail",
                x_le=6.75,
                z_le=0.00,
                num_chordwise_panels=6,
                chordwise_spacing="uniform",
                symmetric=True,
                # Define this wing's root wing cross section.
                wing_cross_sections=[
                    ps.geometry.WingCrossSection(
                        chord=1.5,
                        # Give the root wing cross section an airfoil.
                        airfoil=ps.geometry.Airfoil(
                            name="naca0012",
                        ),
                        twist=-5.0,
                    ),
                    # Define the wing's tip wing cross section.
                    ps.geometry.WingCrossSection(
                        x_le=0.5,
                        y_le=2.0,
                        z_le=1.0,
                        chord=1.0,
                        twist=-5.0,
                        # Give the tip wing cross section an airfoil.
                        airfoil=ps.geometry.Airfoil(
                            name="naca0012",
                        ),
                    ),
                ],
            ),    
        ],
    )
    main_wing_root_wing_cross_section_movement=ps.movement.WingCrossSectionMovement(
        base_wing_cross_section=example_airplane.wings[0].wing_cross_sections[0],
    )
    main_wing_tip_wing_cross_section_movement=ps.movement.WingCrossSectionMovement(
        base_wing_cross_section=example_airplane.wings[0].wing_cross_sections[1],
        sweeping_amplitude=30.0,
        #____________________________________________________________________
        sweeping_period=fp,
        #____________________________________________________________________
        sweeping_spacing="sine",
        pitching_amplitude=0.0,
        pitching_period=0.0,
        pitching_spacing="sine",
        heaving_amplitude=0.0,
        heaving_period=0.0,
        heaving_spacing="sine",
    )
    v_tail_root_wing_cross_section_movement=ps.movement.WingCrossSectionMovement(
        base_wing_cross_section=example_airplane.wings[1].wing_cross_sections[0],
    )
    v_tail_tip_wing_cross_section_movement=ps.movement.WingCrossSectionMovement(
        base_wing_cross_section=example_airplane.wings[1].wing_cross_sections[1],
    )
    main_wing_movement=ps.movement.WingMovement(
        base_wing=example_airplane.wings[0],
        wing_cross_sections_movements=[
            main_wing_root_wing_cross_section_movement,
            main_wing_tip_wing_cross_section_movement,
        ],
    )
    del main_wing_root_wing_cross_section_movement
    del main_wing_tip_wing_cross_section_movement
    v_tail_movement=ps.movement.WingMovement(
        base_wing=example_airplane.wings[1],
        wing_cross_sections_movements=[
            v_tail_root_wing_cross_section_movement,
            v_tail_tip_wing_cross_section_movement,
        ],
    )
    del v_tail_root_wing_cross_section_movement
    del v_tail_tip_wing_cross_section_movement
    airplane_movement=ps.movement.AirplaneMovement(
        base_airplane=example_airplane,
        wing_movements=[main_wing_movement,v_tail_movement],
    )
    del main_wing_movement
    del v_tail_movement
    example_operating_point=ps.operating_point.OperatingPoint(
        #____________________________________________________________________
        density=1.225,
        beta=0.0,
        velocity=va,
        alpha=aoa,
        nu=15.06e-6,
        external_thrust=0.0,
        #____________________________________________________________________
    )
    operating_point_movement=ps.movement.OperatingPointMovement(
        base_operating_point=example_operating_point,
        velocity_amplitude=0.0,
        velocity_period=0.0,
        velocity_spacing="sine",
    )
    movement=ps.movement.Movement(
        airplane_movements=[airplane_movement],
        operating_point_movement=operating_point_movement,
        num_steps=None,
        delta_time=None,
    )
    del airplane_movement
    del operating_point_movement
    example_problem=ps.problems.UnsteadyProblem(
        movement=movement,
    )
    example_solver=ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
        unsteady_problem=example_problem,
    )
    del example_problem
    example_solver.run(
        logging_level="Warning",
        prescribed_wake=True,
    )

    lift, induced_drag, side_force, pitching_moment = extract_forces(example_solver)  

    return lift, induced_drag, pitching_moment


def Mark4Simulation(FlappingPeriod, Angle_of_Attack, Air_Speed):

    FlappingPeriod = [ FlappingPeriod]
    Angle_of_Attack = [Angle_of_Attack]
    Air_Speed = [Air_Speed]
    # Function call
    for i in FlappingPeriod:
        for j in Air_Speed:
            for k in Angle_of_Attack:
                #D Increament Parameter to count the number of CSV Files
                # Call simulation function with Airspeed (k), AoA (j), and Flapping Frequency (i)
                lift, induced_drag, pitching_moment = simulation(j, k, i)

                # Flatten the shit out of them
                lift = lift.flatten()
                induced_drag = induced_drag.flatten()
                pitching_moment = pitching_moment.flatten()
                print("Before Second Cycle Extraction :  ", lift.shape, induced_drag.shape, pitching_moment.shape)
                print("Before Second Cycle Extraction value:  ", lift, induced_drag, pitching_moment)
                
                # Extract the second cycle from the data
                lift, induced_drag, pitching_moment = extract_second_cycle(lift, induced_drag, pitching_moment)
                print("After Second Cycle Extraction :  ",lift.shape, induced_drag.shape, pitching_moment.shape)
                print("After Second Cycle Extraction value:  ", lift, induced_drag, pitching_moment)
                # Save data to CSV with a unique filename for each combination
                file_name = "ModelVerification.csv"
                Mk4SaveDataToCSV(i, j, k, lift, pitching_moment, induced_drag, file_name)

    print("Data Collection is Over")







