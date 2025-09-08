The files above are all related to the results section.

Collision_Avoidance_Maneuver: This file takes the states from the Updating_GS-UKF files and calculates the probabilities of collision, as well as the burst velocity required to reduce probability.

Comparing_GS-EKF_and_Monte_Carlo: This file compares the predictions from GS-EKF and Monte_Carlo predictions, giving the qualitative results.

Comparing_GS-UKF_and_Monte_Carlo: This file compares the predictions from GS-UKF and Monte_Carlo predictions, giving the qualitative results.

Monte_Carlo_ECI: This file gives the propagated distribution using Monte_Carlo, initiated in ECI coordinates.

Monte_Carlo_Equinoctial: This file gives the propagated distribution using Monte_Carlo, initiated in Equinoctial coordinates.

Monte_Carlo_Keplerian: This file gives the propagated distribution using Monte_Carlo, initiated in Keplerian coordinates.

Updating_GS-UKF_Propagated_from_file: This file updates the state propagated with GS-UKF, this technique reads directly from the file.

Updating_GS-UKF_from_elements: This file updates the state propagated with GS-UKF, this technique simply takes the keplerian elements.

Debris_Folder_Total: Move the debris files from this to Debris_Total.

Typically, the ISS.xml has been used as the primary satellite for comparison (note the OMM data isn't actually for the ISS, its just so the name is noticeable).
