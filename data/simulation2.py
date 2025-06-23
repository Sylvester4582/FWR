
import pterasoftware as ps
from SelfFunctions import extract_forces, plot_forces

def sim(va=5.0, aoa=5.0, fp=1):
    fwuav=ps.geometry.Airplane(
        name="FWUAV",
        x_ref=0.0,
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
                        y_le= 0.7,
                        z_le=0.0,
                        chord=0.15,
                        twist=0.0,
                        airfoil=ps.geometry.Airfoil(
                            name="naca8304",
                        ),
                    ),
                ],
            ),
        
            ps.geometry.Wing(
                name="V-Tail",
                x_le=0.52,
                z_le=0.00,
                num_chordwise_panels=6,
                chordwise_spacing="uniform",
                symmetric=True,
            
                wing_cross_sections=[
                    ps.geometry.WingCrossSection(
                        chord=0.23,
                    
                        airfoil=ps.geometry.Airfoil(
                            name="naca0004",
                        ),
                        twist=0.0,
                    ),
                
                    ps.geometry.WingCrossSection(
                        x_le=0.14,
                        y_le=0.17,
                        z_le=0.0,
                        chord=0.03,
                        twist=0.0,
                    
                        airfoil=ps.geometry.Airfoil(
                            name="naca0004",
                        ),
                    ),
                ],
            ),
        ],
    )



    main_wing_root_wing_cross_section_movement = ps.movement.WingCrossSectionMovement(

        base_wing_cross_section=fwuav.wings[0].wing_cross_sections[0],
        sweeping_amplitude=0.0,
        sweeping_period=0.0,
        sweeping_spacing="sine",
        pitching_amplitude=0.0,
        pitching_period=0.0,
        pitching_spacing="sine",
        heaving_amplitude=0.0,
        heaving_period=0.0,
        heaving_spacing="sine",
    )

    main_wing_tip_wing_cross_section_movement = ps.movement.WingCrossSectionMovement(
        base_wing_cross_section=fwuav.wings[0].wing_cross_sections[1],
        sweeping_amplitude=30.0,
        sweeping_period=fp,
        sweeping_spacing="sine",
        pitching_amplitude=15.0,
        pitching_period=0.0,
        pitching_spacing="sine",
        heaving_amplitude=0.0,
        heaving_period=0.0,
        heaving_spacing="sine",
    )

    v_tail_root_wing_cross_section_movement = ps.movement.WingCrossSectionMovement(
        base_wing_cross_section=fwuav.wings[1].wing_cross_sections[0],
    )

    v_tail_tip_wing_cross_section_movement = ps.movement.WingCrossSectionMovement(
        base_wing_cross_section=fwuav.wings[1].wing_cross_sections[1],
    )

    main_wing_movement = ps.movement.WingMovement( 
        base_wing=fwuav.wings[0],

        wing_cross_sections_movements=[
            main_wing_root_wing_cross_section_movement,
            main_wing_tip_wing_cross_section_movement,
        ],

        x_le_amplitude=0.0,
        x_le_period=0.0,
        x_le_spacing="sine",
        y_le_amplitude=0.0,
        y_le_period=0.0,
        y_le_spacing="sine",
        z_le_amplitude=0.0,
        z_le_period=0.0,
        z_le_spacing="sine",
    )

    del main_wing_root_wing_cross_section_movement
    del main_wing_tip_wing_cross_section_movement

    v_tail_movement = ps.movement.WingMovement(

        base_wing=fwuav.wings[1],
        wing_cross_sections_movements=[
            v_tail_root_wing_cross_section_movement,
            v_tail_tip_wing_cross_section_movement,
        ],
    )

    del v_tail_root_wing_cross_section_movement
    del v_tail_tip_wing_cross_section_movement

    airplane_movement = ps.movement.AirplaneMovement(
        base_airplane=fwuav,
        wing_movements=[main_wing_movement, v_tail_movement],
        x_ref_amplitude=0.0,
        x_ref_period=0.0,
        x_ref_spacing="sine",
        y_ref_amplitude=0.0,
        y_ref_period=0.0,
        y_ref_spacing="sine",
        z_ref_amplitude=0.0,
        z_ref_period=0.0,
        z_ref_spacing="sine",
    )

    del main_wing_movement
    del v_tail_movement

    example_operating_point = ps.operating_point.OperatingPoint(
    # -------------------------------------------------------------------------------------------------------------------
        density=1.225,
        beta=0.0,
        velocity=va,
        alpha=aoa,
        nu=15.06e-6,
    # -------------------------------------------------------------------------------------------------------------------
    )
    operating_point_movement = ps.movement.OperatingPointMovement(
        base_operating_point=example_operating_point,
        velocity_amplitude=0.0,
        velocity_period=0.0,
        velocity_spacing="sine",
    )

    movement = ps.movement.Movement(
        airplane_movements=[airplane_movement],
        operating_point_movement=operating_point_movement,
        num_steps=None,
        delta_time=None,
    )

    del airplane_movement
    del operating_point_movement

    example_problem = ps.problems.UnsteadyProblem(
        movement=movement,
    )

    example_solver = ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
        unsteady_problem=example_problem,
    )

    del example_problem

    example_solver.run(
        logging_level="Warning",
        prescribed_wake=True,
    )

    # Animate output

    # ps.output.animate(
    #     unsteady_solver=example_solver,
    #     scalar_type="lift",
    #     show_wake_vortices=True,
    #     save=True,
    # )

    # Extract forces from the solver 

    lift, induced_drag, side_force, pitching_moment = extract_forces(example_solver)
    plot_forces(lift, induced_drag, side_force, pitching_moment)

    return lift, induced_drag, side_force, pitching_moment


sim(5.0, 5.0, 1.0)