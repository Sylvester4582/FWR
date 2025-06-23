
# Importing the header files
import pterasoftware as ps
import numpy as np
from SelfFunctions import print_unsteady_results


# Function to run the unsteady aerodynamic model with conceptual design models and (va, aoa, and flapping period) as the input parameters
def simulation(fp = 0.2, va = 5, aoa = 5):
    
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
                        num_spanwise_panels=6,
                        spanwise_spacing="cosine",
                        chord=0.43,
                        airfoil=ps.geometry.Airfoil(
                            name= "naca0004",
                            coordinates=None,
                            repanel=True,
                            n_points_per_side=400,
                        ),
                    ),
                    ps.geometry.WingCrossSection(
                        x_le=0.0,
                        y_le= 1.4,
                        z_le=0.0,
                        chord=0.15,
                        twist=0.0,
                        airfoil=ps.geometry.Airfoil(
                            name="naca8304",
                        ),
                    ),
                ],
            ),
            # Define the next wing.
            ps.geometry.Wing(
                name="V-Tail",
                x_le=0.43,
                z_le=0.00,
                num_chordwise_panels=6,
                chordwise_spacing="uniform",
                symmetric=True,
                # Define this wing's root wing cross section.
                wing_cross_sections=[
                    ps.geometry.WingCrossSection(
                        chord=0.24,
                        # Give the root wing cross section an airfoil.
                        airfoil=ps.geometry.Airfoil(
                            name="naca0004",
                        ),
                        twist=0.0,
                    ),
                    # Define the wing's tip wing cross section.
                    ps.geometry.WingCrossSection(
                        x_le=0.14,
                        y_le=0.17,
                        z_le=0.0,
                        chord=0.03,
                        twist=0.0,
                        # Give the tip wing cross section an airfoil.
                        airfoil=ps.geometry.Airfoil(
                            name="naca0004",
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
    
    main_wing_movement=ps.movement.WingMovement(
        base_wing=example_airplane.wings[0],
        wing_cross_sections_movements=[
            main_wing_root_wing_cross_section_movement,
            main_wing_tip_wing_cross_section_movement,
        ],
    )
    del main_wing_root_wing_cross_section_movement
    del main_wing_tip_wing_cross_section_movement


    airplane_movement=ps.movement.AirplaneMovement(
        base_airplane=example_airplane,
        wing_movements=[main_wing_movement],
    )

    del main_wing_movement

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

    ps.output.animate(
    # Set the unsteady solver to the one we just ran.
    unsteady_solver=example_solver,
    # Tell the animate function to color the aircraft's wing panels with the local
    # lift coefficient. The valid arguments for this parameter are None, "induced drag",
    # "side force", or "lift".
    scalar_type="lift",
    # Tell the animate function to show the wake vortices. This value defaults to
    # False.
    show_wake_vortices=True,
    # Tell the animate function to not save the animation as file. This way,
    # the animation will still be displayed but not saved. This value defaults to
    # False.
    save=False,
    )

    lift, induced_drag = print_unsteady_results(example_solver)  

    return lift, induced_drag



print(simulation(fp = 0.8, va = 5, aoa = 35))





