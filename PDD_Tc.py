import numpy as np
from time import time


def get_tc(tb, lower_tc, upper_tc, iterations, measured_E, M, variation):

    #   Create random values for tc
    tc = np.random.uniform(lower_tc, upper_tc, iterations)

    #   Calculate erosion rates
    sim_E = (-0.144*(tb/tc)**3+0.904*(tb/tc)**2-0.823*(tb/tc)+0.204)*M*tc

    #   First filter--> remove Tc which tb/tc is lower than 0.52
    tc[tb/tc < 0.52] = 0

    #   Second filter--> remove Tc which shows values higher than the variation selected by the user
    #   variation = simulated erosion rates - measured erosion rates
    tc[abs(sim_E - measured_E) > variation] = 0

    #   Any error just ignore
    np.seterr(invalid='ignore')
    return tc


def read_Photosedfiles(data_txt, time_window, dry_density):
    data_list = []

    #   Open the data text
    data = open(data_txt, 'r')

    for line in data.readlines():
        #   Converts the line into a list using a space separator
        line_as_list = line.split(" ")

        #   Append an empty sub-list for every file line
        data_list.append([])

        for entry in line_as_list:
            try:
                #   Try to append the entry as floating point number to the last sub-list
                #   The last sub-list is pointed at using [-1]

                data_list[-1].append((float(entry)/time_window)*dry_density)
            except ValueError:

                #   If entry is not numeric, append 0.0 to the sub-list and print a warning message
                print("Warning: %s is not a number. Replacing value with 0.0." % str(entry))

    #   Verify that data_list contains the number of rows (sub-lists) with the built-in list function __len__()
    print("Number of rows: %d" % data_list.__len__())

    #   Verify that the first sub-list has the number of columns
    print("Number of columns: %d" % data_list[0].__len__())

    #   Close file (otherwise it will be locked as long as Python is still running!) alternative: use with-statement
    data.close()
    E_rates_array = np.array(data_list)
    num_rows = np.shape(E_rates_array)[0]
    num_cols = np.shape(E_rates_array)[1]

    #   Print the number of rows and columns
    print(num_rows, num_cols)

    #   Print the data
    print(E_rates_array)

    return E_rates_array


def erosion_regions100(array):

    new_array = []
    #   Subdivide the array in small regions upon the user's selection and find the mean value
    #   In this example, the region of interested (ROI) was divided in 49 cells

    mean_matrix1 = np.sum(array[0:100, 0:100])/10000
    mean_matrix2 = np.sum(array[100:200, 0:100])/10000
    mean_matrix3 = np.sum(array[200:300, 0:100])/10000
    mean_matrix4 = np.sum(array[300:400, 0:100])/10000
    mean_matrix5 = np.sum(array[400:500, 0:100])/10000
    mean_matrix6 = np.sum(array[500:600, 0:100])/10000
    mean_matrix7 = np.sum(array[600:700, 0:100])/10000
    mean_matrix8 = np.sum(array[0:100, 100:200])/10000
    mean_matrix9 = np.sum(array[100:200, 100:200])/10000
    mean_matrix10 = np.sum(array[200:300, 100:200])/10000
    mean_matrix11 = np.sum(array[300:400, 100:200])/10000
    mean_matrix12 = np.sum(array[400:500, 100:200])/10000
    mean_matrix13 = np.sum(array[500:600, 100:200])/10000
    mean_matrix14 = np.sum(array[600:700, 100:200])/10000
    mean_matrix15 = np.sum(array[0:100, 200:300])/10000
    mean_matrix16 = np.sum(array[100:200, 200:300])/10000
    mean_matrix17 = np.sum(array[200:300, 200:300])/10000
    mean_matrix18 = np.sum(array[300:400, 200:300])/10000
    mean_matrix19 = np.sum(array[400:500, 200:300])/10000
    mean_matrix20 = np.sum(array[500:600, 200:300])/10000
    mean_matrix21 = np.sum(array[600:700, 200:300])/10000
    mean_matrix22 = np.sum(array[0:100, 300:400])/10000
    mean_matrix23 = np.sum(array[100:200, 300:400])/10000
    mean_matrix24 = np.sum(array[200:300, 300:400])/10000
    mean_matrix25 = np.sum(array[300:400, 300:400])/10000
    mean_matrix26 = np.sum(array[400:500, 300:400])/10000
    mean_matrix27 = np.sum(array[500:600, 300:400])/10000
    mean_matrix28 = np.sum(array[600:700, 300:400])/10000
    mean_matrix29 = np.sum(array[0:100, 400:500])/10000
    mean_matrix30 = np.sum(array[100:200, 400:500])/10000
    mean_matrix31 = np.sum(array[200:300, 400:500])/10000
    mean_matrix32 = np.sum(array[300:400, 400:500])/10000
    mean_matrix33 = np.sum(array[400:500, 400:500])/10000
    mean_matrix34 = np.sum(array[500:600, 400:500])/10000
    mean_matrix35 = np.sum(array[600:700, 400:500])/10000
    mean_matrix36 = np.sum(array[0:100, 500:600])/10000
    mean_matrix37 = np.sum(array[100:200, 500:600])/10000
    mean_matrix38 = np.sum(array[200:300, 500:600])/10000
    mean_matrix39 = np.sum(array[300:400, 500:600])/10000
    mean_matrix40 = np.sum(array[400:500, 500:600])/10000
    mean_matrix41 = np.sum(array[500:600, 500:600])/10000
    mean_matrix42 = np.sum(array[600:700, 500:600])/10000
    mean_matrix43 = np.sum(array[0:100, 600:700])/10000
    mean_matrix44 = np.sum(array[100:200, 600:700])/10000
    mean_matrix45 = np.sum(array[200:300, 600:700])/10000
    mean_matrix46 = np.sum(array[300:400, 600:700])/10000
    mean_matrix47 = np.sum(array[400:500, 600:700])/10000
    mean_matrix48 = np.sum(array[500:600, 600:700])/10000
    mean_matrix49 = np.sum(array[600:700, 600:700])/10000

    #   Construct new array with the mean of each subregions
    new_array = [mean_matrix1, mean_matrix2, mean_matrix3, mean_matrix4, mean_matrix5, mean_matrix6, mean_matrix7, mean_matrix8, mean_matrix9, mean_matrix10, mean_matrix11, mean_matrix12, mean_matrix13, mean_matrix14, mean_matrix15, mean_matrix16, mean_matrix17, mean_matrix18, mean_matrix19, mean_matrix20, mean_matrix21, mean_matrix22, mean_matrix23, mean_matrix24, mean_matrix25, mean_matrix26, mean_matrix27, mean_matrix28, mean_matrix29, mean_matrix30, mean_matrix31, mean_matrix32, mean_matrix33, mean_matrix34, mean_matrix35, mean_matrix36, mean_matrix37, mean_matrix38, mean_matrix39, mean_matrix40, mean_matrix41, mean_matrix42, mean_matrix43, mean_matrix44, mean_matrix45, mean_matrix46, mean_matrix47, mean_matrix48, mean_matrix49]
    E_rates_array = np.array([new_array])

    #   Reshape array to vector
    E_vector = np.reshape(E_rates_array, -1)

    #   Print vector
    print(new_array)
    return E_vector

def main():

    #   Read the files from PHOTOSED--> matrix of eroded depth per cell for specific t (s) at surface erosion mode
    #   The input information are the following: location of the file, time window of the measurement, and dry density of the sediment
    #   In the example, the file used corresponds to homogeneous sediment sample test B (id: AM35) for the last minute at discharge 12 L/s
    e_rates_array = read_Photosedfiles('C:/Users/...../AM35_T06_Q120 543.jpgAM35_T06_Q120 603.jpg_270_sc_bil.txt', 60, 1.32)

    #   Read the mean value of the subregions
    e_vector = erosion_regions100(e_rates_array)

    #   For each value of the subregion, iterate the function get_tc
    for x in np.nditer(e_vector):
        #   The input information for function get_tc are the following: tb, lower_tc, upper_tc, iterations, measured_E, M and variation
        #   tb, lower_tc, upper_tc units are in Pa. E and variation units are in kg/m2*s
        sim_tc = get_tc(1.611732, 0.01, 4, 1000000, x, 0.02, 0.00001)
        start = time()

    #   Filter negative values of tc
        positive_simTc = [x for x in sim_tc if x > 0]
        average_Tc_value = [np.mean(positive_simTc)]

    #   Print number of values that were removed and were used
        print("Points removed", len([i for i in sim_tc if i == 0]))
        print("Points used", 1000000-len([i for i in sim_tc if i == 0]))

        end = time()

    #   Print time of the simulation
        print("time (s)=", end - start)

    #   Create a text file to save average Tc
        average_Tc_txtfile = open('C:/Users/...../SIMTC_AM35_Q12_540-600_Sim10_Me0_02.txt', mode='a')

    #   Write the values in the text file
        average_Tc_txtfile.write(", ".join([str(e) for e in average_Tc_value]) + "\n")

        average_Tc_txtfile.close()


if __name__ == '__main__':
    # launch main function
    main()

