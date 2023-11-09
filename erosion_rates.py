import numpy as np
from time import time


def get_Esim(tb1, tc, M, SD_tb, iterations):

    #   Create random values for tb
    tb = np.random.normal(tb1, SD_tb, iterations)

    #   Create a file to save tb and write them to a text file
    input_file_Tb = open('C:/Users/...../1_Tb_generated.txt', mode='a')
    input_file_Tb.write(str(tb) + "\n")
    input_file_Tb.close()

    #   Calculate simulated erosion rates
    sim_E = (-0.144*(tb/tc)**3+0.904*(tb/tc)**2-0.823*(tb/tc)+0.204)*M*tc

    #   Calculate simulated erosion rates, when tb/tc > 1.7
    sim_E2 = ((tb/tc)-1)*M*tc

    #   First filter--> remove Erosion rates (sim_E) which tb/tc is lower than 0.52 and higher than 1.7
    sim_E[tb / tc < 0.52] = 0
    sim_E[tb / tc > 1.7] = 0

    #   Second filter--> remove Erosion rates (sim_E2) which tb/tc is lower than 1.7
    sim_E2[tb / tc < 1.7] = 0

    #   Sum sim_E and sim_E2. While sim_E2 corresponds to the erosion rates when tb/tc > 1.7, Sim_E is calculated when tb/tc < 1.7 and tb/tc > 0.52
    sim_Etotal = sim_E + sim_E2

    return sim_Etotal


def get_z(sim_Etotal, dt, dryDensity, volfract, zmax):

    #   Calculate the eroded depth (dz)
    dz = (sim_Etotal*dt) / (dryDensity*volfract)
    dz[dz > zmax] = zmax
    return dz


def main():
    #   The following list of variables are selected by the user depending on the equipment and conditions during the measurement
    #   The most important variables are the tc (mean and sd), Q start and Q end
    start = time()
    number_simulations = 1000

    #   Asub, Aroi unit is in m2
    Aroi = 0.002456
    #   The number of iterations is related with the number of subregions. In the example, they were used 49 subregions
    iterations = 49

    #   dry density unit is in kg/m3
    dryDensity = 1320
    volfract = 1

    #   dt refers to the time window for the estimation of the erosion. dt unit is in s
    #   t_end refers to the time window when the bed shear stress (discharge) is kept constant
    dt = 60
    t_end = 600

    #   mean_tc and SD_tc units are in Pa. This values are obtained with the code PDD_Tc
    mean_tc = 1.33
    SD_tc = 0.17
    M = 0.02

    #   zmax, zmaxacum  units are in m. zmax is the maximum eroded depth in time dt
    #   zmaxacum is the maximum eroded depth acumulated over time. In our tests, PHOTOSED was only able to detect depths up to 0.01 m
    zmax = 0.005
    zmaxacum = 0.01

    for i in range(0, number_simulations):

    #   Create random values for tc
        tc = np.random.normal(mean_tc, SD_tc, iterations)

    #   Create dz_acumulate with the number of iterations
        dz_acumulate = np.zeros(iterations)
        vol_subroi = np.zeros(iterations)

    #   Update the information of the zmaxacum
        zmaxacum_subroi = np.zeros(iterations) + zmaxacum

    #   Create a file to save the tc and write them to a text file
        input_file_Tc = open('C:/Users/...../1_Tc_generated.txt', mode ='a')
        input_file_Tc.write(str(tc) + "\n")
        input_file_Tc.close()

    #   Calculate the number of the analyzed time slots for each imposed shear stress (discharge)
        Nt = int(t_end / dt)

    #   Read text files that contain the mean tb and the standard deviation of tb measured at the SETEG channel for each imposed bed shear stress
        tb1_file = np.loadtxt('99_mean_tb_values.txt')
        SD_tb_file = np.loadtxt('99_SD_tb_values.txt')

    #   The iteration is performed for each imposed bed shear stress using the corresponding tb (mean and SD)
        for tb1, SD_tb in zip(tb1_file, SD_tb_file):

            for i in range(0, Nt):
                print("dt->", str(i + 1))
    #   Calculate the erosion rates using the text files of tb and the tc, taking into account the number of iterations
                Esim = np.array(get_Esim(tb1, tc, M, SD_tb, iterations))
    #   Calculate the erosion depth from each simulated erosion rate
                dz = np.array(get_z(Esim, dt, dryDensity, volfract, zmax))
    #   Calculate the accumulated z to compare against the zmax. If the zmax is overpassed, the increment of z is set to 0
                dz_increment = dz
                dz_current = dz_acumulate+dz_increment
                dz_acumulate = dz_current
                values_overZmax = zmaxacum_subroi - dz_acumulate
                values_overZmax[zmaxacum_subroi - dz_acumulate > 0] = 0
                values_overZmax[zmaxacum_subroi - dz_acumulate < 0] = 1
                dz[values_overZmax > 0] = 0
    #   Calculate the volume eroded for each subregion. In the example, they were used 49 subregions
                vol_subroi = dz * (Aroi/iterations)
    #   Calculate the erosion rate for the whole region of interested (ROI)
                ero_rate_ROI = np.sum(vol_subroi)*dryDensity*volfract*(1/(Aroi*dt))
    #   Create a text file to save the dz_acumulate
                filename = str(i+1) + ".txt"
                generated_filename = open('C:/...../dz_acum_t'+filename, mode='a')
                generated_filename.write("\n".join([str(e) for e in dz_acumulate]) + "\n")
                generated_filename.close()

    #   Create a text file to save the simulated erosion rates per cell
                filename = str(i+1) + ".txt"
                generated_filename_sim_E = open('C:/Users/...../sim_ero_cell_t'+filename, mode='a')
                generated_filename_sim_E.write("\n".join([str(e) for e in Esim]) + "\n")
                generated_filename_sim_E.close()

    #   Create a text file to save the simulated erosion rates for the whole ROI
                filename_sim_ero = open('C:/Users/...../sim_erosion_ROI.txt', mode='a')
                filename_sim_ero.write(str(ero_rate_ROI) + "\n")
                filename_sim_ero.close()

    #   Create a text file to save the simulated erosion rates for the whole ROI related with the total duration of the erosion test
                filename_sim_ero = open('C:/Users/...../0_sim_erosion_ROIall.txt', mode='a')
                filename_sim_ero.write(str(ero_rate_ROI) + "\n")
                filename_sim_ero.close()

#   Estimate perfomance indicators of the simulations
#   Comparisson of the simulated erosion rates against measured rates, between Qstart and Qend
        Qstart = 100
        Qend = 110
        num_values = Qend - Qstart

#   Read the file with the measured erosion rates. In the example, the file of the measurements corresponds to an homogeneous mixtures test C (id: AM26)
        measured_E = 'C:/Users/...../99_measured_E_ROI_AM26.txt'
        data_measured_E1 = np.loadtxt(measured_E)
        data_measured_E = data_measured_E1[Qstart:Qstart+num_values]

#   Read the file with the simulated erosion rates
        data_sim_E1 = np.loadtxt('C:/Users/...../sim_erosion_ROI.txt')
        data_sim_E = data_sim_E1[Qstart:Qstart+num_values]
#   Calculate residuals and the mean from the measured erosion rates
        residuals = data_sim_E - data_measured_E
        mean_observed = np.mean(data_measured_E)
#   Calculate performance indicators: RMSE, R-squared, Nash Coefficient, NPE, MAE and KGE
        #   Calculate RMSE
        rmse = np.sqrt(np.mean(residuals ** 2))

        #   Write the results of RMSE to a text file
        output_file1 = open('C:/Users/...../0_RMSE.txt', mode ='a')
        output_file1.write(str(rmse) + "\n")
        output_file1.close()

        #   Calculate R-squared
        ss_total = np.sum((data_measured_E - mean_observed) ** 2)
        ss_residual = np.sum(residuals ** 2)
        r_squared = 1 - (ss_residual / ss_total)

        #   Write the results of R-squared to a text file
        output_file2 = open('C:/Users/...../0_r2.txt', mode ='a')
        output_file2.write(str(r_squared) + "\n")
        output_file2.close()

        #   Calculate Nash Coefficient
        nash = 1 - (np.sum(residuals ** 2) / np.sum((data_measured_E - mean_observed) ** 2))

        #   Write the results of Nash Coefficient to a text file
        output_file3 = open('C:/Users/...../0_nash.txt', mode ='a')
        output_file3.write(str(nash) + "\n")
        output_file3.close()

        #   Calculate Normalized Peak Error (NPE)
        npe = np.max(np.abs(residuals)) / (np.max(data_measured_E) - np.min(data_measured_E))

        #   Write the results of NPE to a text file
        output_file4 = open('C:/Users/...../0_npe.txt', mode ='a')
        output_file4.write(str(npe) + "\n")
        output_file4.close()

        #   Calculate Mean Absolute Error (MAE)
        mae = np.mean(np.abs(residuals))

        #   Write the results of MAE to a text file
        output_file5 = open('C:/Users/...../0_mae.txt', mode ='a')
        output_file5.write(str(mae) + "\n")
        output_file5.close()

        #   Calculate Kling Gupta Efficiency (KGE)
        kge = 1 - np.sqrt((np.corrcoef(data_measured_E, data_sim_E)[0, 1] - 1) ** 2 +
                          (np.std(data_sim_E) / np.std(data_measured_E) - 1) ** 2 +
                          (np.mean(data_sim_E) / np.mean(data_measured_E) - 1) ** 2)

        # Write the results of KGE to a text file
        output_file6 = open('C:/Users/...../0_kge.txt', mode ='a')
        output_file6.write(str(kge) + "\n")
        output_file6.close()

        #   Calculate AVERAGE values
        #   Write the average Erosion rate from one simulation
        output_file7 = open('C:/Users/...../0_average_simE.txt', mode ='a')
        output_file7.write(str(np.average(data_sim_E)) + "\n")

        measured_E_averaged = open('C:/Users/...../99_measured_E_ROI_ave.txt', mode ='a')
        measured_E_averaged.write(str(np.average(data_measured_E)) + "\n")

        #   Write down the file path of the simulated erosion rates
        file_path ='C:/Users/...../sim_erosion_ROI.txt'
        with open(file_path,'w'):
            pass

    end = time()
    print("time (s)=", end - start)

if __name__ == '__main__':
    # launch main function
    main()