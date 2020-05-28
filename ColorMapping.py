import sys
import math
import collections
from matplotlib import cm
import matplotlib.colors as colors
import matplotlib.pylab as pltlab
import numpy
import gdal

dem_name = "/Users/alex/Downloads/HerodotusA_A_dem.cub"


def color_mapping(dem_name, a, *args):
    dem_file = gdal.Open(dem_name)
    dem = numpy.array(dem_file.GetRasterBand(1).ReadAsArray())
    dem[dem == -3.40282265508890445e+38] = numpy.NaN
    log_it = open(dem_name + ' log.txt', 'w')
    if a < 0:
        print('\nERROR: a = alpha value. It must be between 0 (inclusive) and 90 (exclusive)\n')
        log_it.write('\nERROR: a = alpha value. It must be between 0 (inclusive) and 90 (exclusive)\n')
        sys.exit()
    elif a > 90:
        print('\nERROR: a = alpha value. It must be between 0 (inclusive) and 90 (exclusive)\n')
        log_it.write('\nERROR: a = alpha value. It must be between 0 (inclusive) and 90 (exclusive)\n')
        sys.exit()
    elif a == 90:
        a = 89.99999999999999
    # alpha values:
    # a=89.9999999999; %angle between 0 and 90.
    # a=45
    # a=1
    # a=45; %input alpha value

    # Import Color Map

    if 'lut' in str.lower(str(args)):
        print('\nImporting LUT file\n')
        log_it.write('\nImporting LUT file\n')
        # TODO Read in LUT file
    elif ' m ' in str.lower(str(args)):
        print('Importing matlab colormap')
        log_it.write('\nImporting matlab colormap\n')
        # TODO assign map from colormap file
    try:
        color_map_type = cm.get_cmap(args[1])
        norm_map_alpha = numpy.array(color_map_type.colors)
        print('Using ' + args[1] + ' Colormap')
        log_it.write('Using ' + args[1] + ' Colormap')
        norm_map = norm_map_alpha[:, :3]
        map = norm_map * 255
        print(map)
    except ValueError:
        print('No lut, matlab colormap, or regular colormap; using standard map')
        log_it.write('No lut or matlab colormap using standard map')
        map = numpy.array([[30, 83, 151],
                        [69, 117, 180],
                        [97, 144, 195],
                        [125, 172, 209],
                        [154, 197, 223],
                        [185, 217, 234],
                        [215, 238, 245],
                        [233, 247, 233],
                        [255, 255, 191],
                        [255, 244, 174],
                        [255, 233, 157],
                        [254, 215, 138],
                        [253, 183, 116],
                        [253, 150, 95],
                        [243, 116, 76],
                        [229, 82, 57],
                        [215, 48, 39]])
        norm_map = map / 255
    except IndexError:
        print('No lut, matlab colormap, or regular colormap; using standard map')
        log_it.write('No lut or matlab colormap using standard map')
        map = numpy.array([[30, 83, 151],
                           [69, 117, 180],
                           [97, 144, 195],
                           [125, 172, 209],
                           [154, 197, 223],
                           [185, 217, 234],
                           [215, 238, 245],
                           [233, 247, 233],
                           [255, 255, 191],
                           [255, 244, 174],
                           [255, 233, 157],
                           [254, 215, 138],
                           [253, 183, 116],
                           [253, 150, 95],
                           [243, 116, 76],
                           [229, 82, 57],
                           [215, 48, 39]])
        norm_map = map / 255

    # get dimensions and sort values
    # ordered values
    log_it.write(str(dem.shape))
    [m, n] = numpy.size(dem, 0), numpy.size(dem, 1)
    log_it.write('\n[m, n]=' + str([m, n]) + '\n')
    D = numpy.around(numpy.sort(dem, axis=None), decimals=10)
    I = numpy.argsort(dem, axis=None)
    print(I[0])
    # log_it.write('\nD= ' + str(D) + '\n')
    # log_it.write('\nI= ' + str(I) + '\n')

    # Finding max and min Index non null number
    max_elev = numpy.nanmax(D)
    min_elev = numpy.nanmin(D)
    max_idx = (numpy.count_nonzero(~numpy.isnan(D)))-1

    # EQUATION #
    # Let v_diag = v_n - v_1 and p_i = < v_i, v_diag>/|v_diag|^2
    # Data Driven Color Mapping by Martin Eisemann, Georgia Albuquerque, Marcus
    # Magnor
    # https://pdfs.semanticscholar.org/1145/6c8c21cd2b740561b6b82516fd61ce125051.pdf

    # x is position in sorted array
    # y is data value

    # To change angle of diagonal:
    # vdiag = vn - v1 = [xn,yn] - [x1,y1] = [xn-x1,yn-y1] = [dx,dy]
    # dx = (yn-y1)/tan(a)

    # Data Projection
    # Calculate v_diag & V:
    # dx = max_idx - 1
    try:
        dx = numpy.array((D[max_idx] - D[0]) / math.tan(math.radians(a)))
    except ZeroDivisionError:
        print("encountered zero division error error in determining dx")
        log_it.write("encountered zero division error in determining dx")
        sys.exit()
    dx_ = numpy.repeat(dx, len(D))
    dy = numpy.array(D[max_idx] - D[0])
    dy_ = numpy.repeat(dy, len(D))
    v_diag = [dx, dy]

    # Square of magnitude of v_diag.
    denominator = numpy.linalg.norm(v_diag) * numpy.linalg.norm(v_diag)
    
    v_ida_numerator = numpy.array(range(0, len(D)))
    v_idx = numpy.array(v_ida_numerator / len(D) * dx)
    #  print('\nmin v_idx = ' + str(numpy.nanmin(v_idx)) + '\n')
    #  print('\nmax v_idx = ' + str(numpy.nanmax(v_idx)) + '\n')
    # vector of positions (indices) (why??)
    # v_idx=( 1:1:length(D) ) ; %vector of positions (indices)
    v_i = numpy.column_stack((v_idx, D))  # 2 by n matrix v_i

    V_diag = numpy.column_stack((dx_, dy_))

    # P = numpy.array(numpy.einsum("ij,ij->i", v_i, V_diag) / denominator)
    P = (numpy.sum(v_i*V_diag, axis=1) / denominator)
    # min(P)
    # max(P)
    print('min P = ' + str(numpy.nanmin(P)) + '\n')
    print('max P = ' + str(numpy.nanmax(P)) + '\n')
    log_it.write('min P = ' + str(numpy.nanmin(P)) + '\n')
    log_it.write('max P = ' + str(numpy.nanmax(P)) + '\n')

    # Put values back into Matrix
    C = [x for x in I]
    P_pmin = numpy.array(P - numpy.nanmin(P))
    Imv = int(numpy.nanmin(I))
    while Imv < numpy.size(I):
        I_index_value = I[Imv]
        C[I_index_value] = P_pmin[Imv]
        Imv += 1
    color_dem = numpy.reshape(C, (m, n))

    # Colormaps

    # Starry Night
    # sn_map=[0,0,102;60,83,163;94,192,238;174,245,255;255,255,255;253,245,124;251,217,40;240,165,0;200,90,0;162,40,10;88,1,1];
    # sn_norm_map=sn_map/255;
    # 
    # Eyeball
    # eb_map=[109,61,103;218,122,206;190,110,158;170,113,217;134,79,182;81,30,145;32,16,200;16,64,200;16,110,220;32,150,240;32,190,255;32,224,255;90,225,245;60,190,220;120,190,220;170,190,220;210,220,210;150,220,150;66,200,75;90,175,75;120,150,85;100,100,60;130,80,20;160,100,25;190,130,30;230,170,20;255,239,10;250,140,0;230,70,0;255,0,0;250,95,95;250,180,155;250,250,250];
    # eb_norm_map=eb_map/255;

    # ColorBrewer
    # map=[30,83,151;69,117,180;97,144,195;125,172,209;154,197,223;185,217,234;215,238,245;233,247,233;255,255,191;255,244,174;255,233,157;254,215,138;253,183,116;253,150,95;243,116,76;229,82,57;215,48,39];
    # norm_map=map/255;
    # 
    # %Normal
    # n_map=[255,120,255;120,120,255;120,255,255;120,255,120;255,255,120;255,120,120;255,255,255];
    # n_norm_map=n_map/255;

    # Viridis (17) norm_map=[0.267004, 0.004874, 0.329415;0.282327, 0.094955, 0.417331;0.278826, 0.175490,
    # 0.483397;0.258965, 0.251537, 0.524736;0.229739, 0.322361, 0.545706;0.199430, 0.387607, 0.554642;0.172719,
    # 0.448791, 0.557885;0.149039, 0.508051, 0.557250;0.127568, 0.566949, 0.550556;0.120638, 0.625828,
    # 0.533488;0.157851, 0.683765, 0.501686;0.246070, 0.738910, 0.452024;0.369214, 0.788888, 0.382914;0.515992,
    # 0.831158, 0.294279;0.678489, 0.863742, 0.189503;0.845561, 0.887322, 0.099702;0.993248, 0.906157, 0.143936];
    # map=norm_map*255;

    # %Magma (18)0.001462, 0.000466, 0.013866; norm_map=[0.035520, 0.028397, 0.125209;0.102815, 0.063010,
    # 0.257854;0.191460, 0.064818, 0.396152;0.291366, 0.064553, 0.475462;0.384299, 0.097855, 0.501002;0.475780,
    # 0.134577, 0.507921;0.569172, 0.167454, 0.504105;0.664915, 0.198075, 0.488836;0.761077, 0.231214,
    # 0.460162;0.852126, 0.276106, 0.418573;0.925937, 0.346844, 0.374959;0.969680, 0.446936, 0.360311;0.989363,
    # 0.557873, 0.391671;0.996580, 0.668256, 0.456192;0.996727, 0.776795, 0.541039;0.992440, 0.884330,
    # 0.640099;0.987053, 0.991438, 0.749504]; map=norm_map*255; Morgan's

    # Display Images
    # TODO probably move this below final? maybe option to display both?

    if 'display' in args or 'Display' in args:
        color_map = colors.ListedColormap(norm_map)
        pltlab.imshow(color_dem)
        pltlab.set_cmap(color_map)
        pltlab.show()

    # Initialize Elevation/Bin related Variables

    elev_range = max_elev - min_elev
    step = (numpy.nanmax(P_pmin) - numpy.nanmin(P_pmin)) / (numpy.size(map, axis=0))
    elevs = numpy.array(numpy.zeros((int(numpy.size(map, axis=0)))))

    # Figure out preliminary elevation values
    # add min and max elev and pcts
    elevs[0] = min_elev
    # fill in the rest of the elevations
    l = 1
    while l < numpy.size(map, axis=0):  # 17
        #  P_pmin = numpy.around(numpy.array(P - numpy.nanmin(P)), decimals=15)
        step_l = numpy.around((step * l), decimals=15)
        P_pmin_step_l = (P_pmin - step_l)
        absolute_P_step = numpy.around(numpy.absolute(P_pmin_step_l), decimals=15)
        Min_absolute_P_step = numpy.nanmin(absolute_P_step)
        closestIndex = numpy.array(numpy.where(absolute_P_step == Min_absolute_P_step))
        elevs[l] = D[closestIndex.flat[-1]]
        l += 1

    # Figure out elevation Differences, and Min Bin Size
    elev_diffs = numpy.array(numpy.diff(elevs, n=1, axis=0))
    min_diff = numpy.nanmin(elev_diffs)  # TODO compare this to min bin size for final rounding
    max_diff = numpy.nanmax(elev_diffs)
    # rounding options
    round_opts = numpy.array([1, 5, 10, 25, 50, 100])  # potential ways to round can be edited
    # try and find out which minimum bin size is closet to the size of the smallest difference
    round_diffs = round_opts - min_diff
    round_diffs = [math.inf if x <= -1 else x for x in round_diffs]
    round_I = list(numpy.where(numpy.absolute(round_diffs) == numpy.nanmin(
        numpy.absolute(round_diffs))))  # rounding option that's > or ~= to minimum

    # Make sure min_bin_size isn't larger than max_bin_size
    while max_diff < round_opts[tuple(round_I)]:
        round_I -= round_I

    # Make sure legend can handle minimum bin size
    min_percent = 0.02 * elev_range
    if min_percent > round_opts[tuple(round_I)]:
        min_percent_diffs = numpy.ndarray.flatten(round_opts) - min_percent
        min_percent_diffs = [math.inf if x <= -1 else x for x in min_percent_diffs]
        round_I = list(numpy.where(
            min_percent_diffs == numpy.nanmin(min_percent_diffs)))  # rounding option that's > or != to minimum
        min_value = min_percent
    else:
        min_value = min_diff

    min_bin_size = 0  # Must define variable outside of statement before changing it
    if round_I[0] == 0:
        # NEAREST 1
        min_bin_size = round(min_value)
        if min_bin_size == 0:
            min_bin_size = 1
    elif round_I[0] == 1:
        # NEAREST 5:
        min_bin_size = round(min_value / 5) * 5
        if min_bin_size == 0:
            min_bin_size = 5
    elif round_I[0] == 2:
        # NEAREST 10:
        min_bin_size = round(min_value, -1)
        if min_bin_size == 0:
            min_bin_size = 10
    elif round_I[0] == 3:
        # NEAREST 25:
        min_bin_size = round(min_value / 25) * 25
        if min_bin_size == 0:
            min_bin_size = 25
    elif round_I[0] == 4:
        # NEAREST 50:
        min_bin_size = round(min_value / 50) * 50
        if min_bin_size == 0:
            min_bin_size = 50
    elif round_I[0] == 5:
        # NEAREST 100:
        min_bin_size = round(min_value, -2)
        if min_bin_size == 0:
            min_bin_size = 100
    else:
        print('Oh no! Something went wrong with determining how to round\n')
        log_it.write('Oh no! Something went wrong with determining how to round\n')
        sys.exit()
    # Which Bins are too small - Initial pass
    S = collections.deque([], maxlen=numpy.size(elev_diffs))  # Stack for center bin numbers
    Dir = collections.deque([], maxlen=numpy.size(elev_diffs))  # stack for direction to shift bins
    # Deque's are a type of stack in python that can be pushed or popped from either side

    # TODO maybe too_small is not necessary
    too_small = (elev_diffs < min_bin_size)  # logical array of bins that are too small 1 = true
    too_small = numpy.multiply(too_small, 1)
    k = 0
    while k < numpy.size(elev_diffs):
        if too_small[k] == 1:  # find an instance of the first in a streak of too small bins
            j = k  # elev_diffs index
            b = 0  # Number of too small bins
            while j < numpy.size(elev_diffs):  # TODO make this an if statement
                if too_small[j] == 0:  # find the last instance of the streak of too small bins
                    # Create a stack (FILO) of the center bins and the
                    # directions that they need to be adjusted:
                    # -1: Adjust the Bin Minimum
                    # 0: Adjust both Bin Min and Max
                    # 1: Adjust the Bin Maximum
                    if b % 2 == 0:  # even number of too small bins
                        center_bin_big = j - b + (b / 2)
                        center_bin_sm = center_bin_big - 1
                        S.append(center_bin_big)
                        Dir.append(1)
                        S.append(center_bin_sm)
                        Dir.append(-1)
                    else:  # odd number of too small bins
                        center_bin_num = j - b + (b // 2)
                        S.append(center_bin_num)
                        Dir.append(0)

                    k = j  # This exits the while loop
                b += 1  # TODO check that this is correct
                j += 1
        k += 1

    # adjust bin sizes
    # Initialize array for adjusted elevs
    elevs_adj = [x for x in elevs]

    while S:
        bin_num = S.pop()
        expandDir = Dir.pop()
        if expandDir == 0:
            # TODO - add check if max bin, so following if doesn't break??
            if (elevs_adj[bin_num + 1] - (elevs_adj[bin_num])) < min_bin_size:
                if bin_num == 1:  # min edge case
                    S.append(bin_num)
                    Dir.append(1)
                elif bin_num == numpy.size(elevs_adj, 2) - 1:
                    S.append(bin_num)
                    Dir.append(-1)
                else:
                    binCenter = ((elevs_adj[bin_num + 1]) + (elevs_adj[bin_num])) / 2
                    new_min = binCenter - (min_bin_size / 2)
                    if new_min <= min_elev:
                        elevs_adj[bin_num] = min_elev
                        S.append(bin_num)
                        Dir.append(1)
                    else:
                        elevs_adj[bin_num] = new_min
                        S.append(bin_num - 1)
                        Dir.append(-1)

                    new_max = (binCenter + (min_bin_size / 2))

                    if new_max >= max_elev:
                        elevs_adj[bin_num + 1] = max_elev
                        S.append(bin_num)
                        Dir.append(-1)
                    else:
                        elevs_adj[bin_num + 1] = new_max
                        S.append(bin_num + 1)
                        Dir.append(1)
        elif expandDir == 1:
            if (elevs_adj[int(bin_num + 1)]) - elevs_adj[int(bin_num)] < min_bin_size:
                if bin_num == float(int(numpy.size(elevs_adj, 0)) - 1):
                    S.append(int(bin_num))
                    Dir.append(-1)
                else:
                    new_max = elevs_adj[int(bin_num)] + min_bin_size
                    if new_max >= max_elev:
                        elevs_adj[int(bin_num + 1)] = max_elev
                        S.append(int(bin_num + 1))
                        Dir.append(-1)
                    else:
                        elevs_adj[int(bin_num + 1)] = new_max
                        S.append(int(bin_num + 1))
                        Dir.append(1)
        elif expandDir == -1:
            if (elevs_adj[int(bin_num + 1)] - elevs_adj[int(bin_num)]) < min_bin_size:
                if bin_num == 1:  # min edge case
                    S.append(int(bin_num - 1))
                    Dir.append(1)
                else:
                    new_min = elevs_adj[int(bin_num+1)] - min_bin_size
                    if new_min <= min_elev:
                        elevs_adj[bin_num] = min_elev
                        S.append(int(bin_num-1))
                        Dir.append(1)
                    else:
                        elevs_adj[int(bin_num)] = new_min
                        S.append(int(bin_num - 1))
                        Dir.append(-1)
        else:
            print("\nOh god. something bad has happened you might what to fix that")
            log_it.write("\nOh god. something bad has happened you might what to fix that\n")

    # elevs_adj /= 1000
    # Round Results

    final_elevs = [x for x in elevs_adj]  # Initial array
    max_round_bin = int(numpy.size(elevs_adj, axis=0) - 1)
    Smallest_final_elevs = final_elevs[0] # not to be rounded so will be replaced after rounding
    Largest_final_elevs = final_elevs[-1]
    print(final_elevs)
    print(max_round_bin)
    if round_I[0] == 0:
        # NEAREST 1
        final_elevs = [round(x) for x in final_elevs]
    elif round_I[0] == 1:
        # NEAREST 5:
        final_elevs = [(x / 5) for x in final_elevs]
        final_elevs = [round(x) for x in final_elevs]
        final_elevs = [(x * 5) for x in final_elevs]
    elif round_I[0] == 2:
        # NEAREST 10:
        final_elevs = [round(x, -1) for x in final_elevs]
    elif round_I[0] == 3:
        # NEAREST 25:
        final_elevs = [(x / 25) for x in final_elevs]
        final_elevs = [round(x) for x in final_elevs]
        final_elevs = [(x * 25) for x in final_elevs]
    elif round_I[0] == 4:
        # NEAREST 50:
        final_elevs = [(x / 50) for x in final_elevs]
        final_elevs = [round(x) for x in final_elevs]
        final_elevs = [(x * 50) for x in final_elevs]
    elif round_I[0] == 5:
        # NEAREST 100:
        final_elevs = [round(x, -2) for x in elevs_adj]
    else:
        print('Oh no! Something went wrong with determining how to round the results')
        log_it.write('\nOh no! Something went wrong with determining how to round the results\n')
    final_elevs = numpy.array(final_elevs)
    # final_elevs /= 1000
    final_elevs[0] = Smallest_final_elevs
    final_elevs[-1] = Largest_final_elevs
    # calculate percentages and make output file
    log_it.write('\nClosing log... \nCreating lut file\n')
    log_it.close()
    print(final_elevs)
    pcts_pre_pre = (final_elevs - min_elev)
    pcts_pre = pcts_pre_pre / elev_range
    pcts = pcts_pre * 100
    fileID = open(dem_name + ' .lut', 'w')
    fileID.write('nv 0 0 0 //noData to black\n')  # NO VALUE LINE
    fileID.write('0% ' + str(map[0, 0]) + ' ' + str(map[0, 1]) + ' ' + str(map[0, 2]) + '\n')  # 0% Line
    for e in range(1, numpy.size(pcts, 0)):
        fileID.write(str(numpy.around(pcts[e], decimals=5)) + '% ' + str(map[(e - 1), 0]) + ' ' + str(map[(e - 1), 1]) + ' ' + str(map[(e - 1), 2]) + '\n')
        fileID.write(str(numpy.around(pcts[e], decimals=5)) + '% ' + str(map[e, 0]) + ' ' + str(map[e, 1]) + ' ' + str(map[e, 2]) + '\n')
    fileID.write('100% ' + str(map[numpy.size(map, 0)-1, 0]) + ' ' + str(map[numpy.size(map, 0)-1, 1]) + ' ' + str(map[numpy.size(map, 0)-1, 2]) + '\n')  # 100%
    fileID.close()
    print('Completed ')
    print(dem_name)
    print('\n')


color_mapping(dem_name, 45, 'Display', 'plasma')
