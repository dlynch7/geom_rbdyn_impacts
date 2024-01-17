# geom_rbdyn_impacts
This repo combines _event detection_ and a geometric integrator for coordinate-free simulation of rigid body dynamics with impacts.

This repo is focused around a single example: a bouncing oblong ball (for example, an American football).
This is a reasonable application for coordinate-free simulation because the ball can tumble and bounce wildly enough that no one coordinate chart is likely to be sufficient for representing the ball's configuration for all time.

This repo uses a 4th-order Crouch-Grossman geometric integrator; see my [other repo](https://github.com/dlynch7/geometric_rigid_body_dynamics) for details on this integrator.

Dependencies: [Modern Robotics](https://github.com/NxRLab/ModernRobotics).