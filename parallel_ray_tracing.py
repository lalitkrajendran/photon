def generate_lightfield_angular_data(lens_pitch, image_distance, scattering_data, scattering_type, lightfield_source,
                                     lightray_number_per_particle, n_min, n_max):
    # % This function generates the lightfield data for the source points specified by the
    # % structure lightfield_source.  The data is only generated for the source points from
    # % n_min to n_max.  The parameter lightray_number_per_particle is the number of rays to generate for each
    # % source point.


    # % If the scattering type is 'mie' then the Mie scattering data is loaded,
    # % otherwise nothing is loaded from the scattering data
    if scattering_type == 'mie':
        #% This saves the scattering angle data into the parameters structure
        mie_scattering_angle = scattering_data['scattering_angle']
        #% This saves the scattering irradiance values for the different particle
        #% diameters into the parameters structure
        mie_scattering_irradiance=scattering_data['scattering_irradiance']
        #% This extracts the inverse rotation matrix from the parameters structure
        inverse_rotation_matrix=scattering_data['inverse_rotation_matrix']
        #% This extracts the normalized beam propogation direction vector from the
        #% parameters structure
        beam_propogation_vector=scattering_data['beam_propogation_vector']

    # % This is the number of points to calculate the lightfield structure for
    source_point_number = n_max - n_min + 1

    # % This initializes the data structure
    lightfield_data = {}
    lightfield_data['x'] = np.zeros(source_point_number * lightray_number_per_particle)
    lightfield_data['y'] = np.zeros(source_point_number * lightray_number_per_particle)
    lightfield_data['z'] = np.zeros(source_point_number * lightray_number_per_particle)
    lightfield_data['theta'] = np.zeros(source_point_number * lightray_number_per_particle)
    lightfield_data['phi'] = np.zeros(source_point_number * lightray_number_per_particle)
    lightfield_data['radiance'] = np.zeros(source_point_number * lightray_number_per_particle)

    # % This is the radius of the current lens for which the light rays are being
    # % generated
    R = lens_pitch / 2.0

    # % This iterates through the "particle" locations generating the light rays
    for n in range(n_min,n_max+1):
        # % This is the current source point for which the lightfield is being generated
        x_current = np.squeeze(np.asarray(lightfield_source['x']))[n]
        y_current = np.squeeze(np.asarray(lightfield_source['y']))[n]
        z_current = np.squeeze(np.asarray(lightfield_source['z']))[n]

        # % This creates random radial coordinates for the lightrays to intersect
        # % on the lens
        np.random.seed(1105)
        r_temp = np.random.rand(1, lightray_number_per_particle)
        r = R * np.sqrt(r_temp)

        # % This creates random angular coordinates for the lightrays to
        # % intersect on the lens
        np.random.seed(4092)
        psi_temp = np.random.rand(1, lightray_number_per_particle)
        psi = 2.0 * np.pi * psi_temp

        # % This calculates the random cartesian coordinate of the points the
        # % rays will intersect on the lens
        x_lens = r * np.cos(psi)
        y_lens = r * np.sin(psi)


        # % This calculates the angular data using Mie scattering if specified in the
        # % parameters structure, otherwise the particles are assummed to have
        # % uniform irradiance
        if scattering_type=='mie':

            # % This extracts the current particle diameter index
            diameter_index=lightfield_source['diameter_index'][n]

            # % This calculates the lightrays direction vectors (in the camera
            # % coordiante system)
            ray_direction_vector = np.array([np.squeeze(x_lens-x_current),np.squeeze(y_lens-y_current),image_distance*np.ones(lightray_number_per_particle)-np.squeeze(z_current)])

            # % This normalizes the ray direction vectors
            ray_direction_vector /=np.sqrt(ray_direction_vector[0,:]**2+ray_direction_vector[1,:]**2+ray_direction_vector[2,:]**2)[np.newaxis,:]

            # % This rotates the lightrays direction vectors by the inverse of the
            # % camera rotation array so that the ray is now in the world coordinate
            # % system
            ray_direction_vector=inverse_rotation_matrix*ray_direction_vector
            # % This calculates the angle that the light ray direction vectors make
            # % with the laser propogation direction vector
            ray_scattering_angles=np.squeeze(np.arccos(beam_propogation_vector*ray_direction_vector))

            # % This calculates the Mie scattering irradiances at the currently
            # % scattered angles and with the current particle diameter
            set_interp = interp1d(mie_scattering_angle,mie_scattering_irradiance[:,diameter_index])
            ray_scattering_irradiance = np.squeeze(set_interp(ray_scattering_angles))
            # ray_scattering_irradiance=interp1d(mie_scattering_angle,mie_scattering_irradiance(:,diameter_index),ray_scattering_angles,'linear')

            # % This calculates the total irradiance for the current particle's rays
            irradiance_current=ray_scattering_irradiance*lightfield_source['radiance'][n]

        elif scattering_type =='diffuse':

            # % This specifies the total irradiance for the current particle's
            # % rays to be uniform
            irradiance_current = lightfield_source['radiance'][n]

        # REMOVE THIS LINE IF YOU REMOVE THE BLOCK COMMENTS
        #irradiance_current = lightfield_source['radiance'][n]

        # % This calculates the x angles for the light rays
        theta_temp = -(np.squeeze(x_lens)- x_current) / (image_distance - z_current)
        # % This calculates the y angles for the light rays
        phi_temp = -(np.squeeze(y_lens) - y_current) / (image_distance - z_current)

        # % This is the vector of indicies to save the light rays to
        index_vector = np.array((n - n_min) * lightray_number_per_particle + np.arange(0, lightray_number_per_particle)).astype('int')
        # index_vector -= 1

        # % This saves the light rays to the data structure
        lightfield_data['x'][index_vector] = x_current
        lightfield_data['y'][index_vector] = y_current
        lightfield_data['z'][index_vector] = z_current
        lightfield_data['theta'][index_vector] = theta_temp
        lightfield_data['phi'][index_vector] = phi_temp
        lightfield_data['radiance'][index_vector] = irradiance_current

    return lightfield_data
