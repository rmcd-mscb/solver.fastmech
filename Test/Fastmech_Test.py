from __future__ import print_function  # Only Python 2.x
import numpy as np
import matplotlib.pyplot as plt
from shutil import copyfile
import h5py
import vtk
import subprocess
import os
from itertools import count
import configparser
import sys

def getCellValue(vtkSGrid2D, newPoint2D, cellID, valarray):
    pcoords = [0.0,0.0,0.0]
    weights = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    clspoint = [0.,0.,0.]
    tmpid = vtk.mutable(0)
    vtkid2 = vtk.mutable(0)
    vtkcell2D = vtk.vtkQuad()
    vtkcell2D = vtkSGrid2D.GetCell(cellID)
    tmpres = vtkcell2D.EvaluatePosition(newPoint2D,clspoint,tmpid,pcoords,vtkid2, weights )
    print(newPoint2D, clspoint, tmpid, pcoords, vtkid2, weights)
    idlist1 = vtk.vtkIdList()
    numpts = vtkcell2D.GetNumberOfPoints()
    idlist1 = vtkcell2D.GetPointIds()
    tmpVal = 0.0
    for x in range(0,numpts):
        tmpVal = tmpVal + weights[x]*valarray.GetTuple(idlist1.GetId(x))[0]
    return tmpVal

def isCellWet(vtkSGrid2D, newPoint2D, cellID, IBC_2D):
    pcoords = [0.0,0.0,0.0]
    weights = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    clspoint = [0.,0.,0.]
    tmpid = vtk.mutable(0)
    vtkid2 = vtk.mutable(0)
    vtkcell2D = vtk.vtkQuad()
    vtkcell2D = vtkSGrid2D.GetCell(cellID)
    tmpres = vtkcell2D.EvaluatePosition(newPoint2D,clspoint,tmpid,pcoords,vtkid2, weights )
    idlist1 = vtk.vtkIdList()
    numpts = vtkcell2D.GetNumberOfPoints()
    idlist1 = vtkcell2D.GetPointIds()
    tmpIBC = 0.0
    for x in range(0,numpts):
        tmpIBC = tmpIBC + weights[x]*abs(IBC_2D.GetTuple(idlist1.GetId(x))[0])
        # print(tmpIBC,abs(IBC_2D.GetTuple(idlist1.GetId(x))[0]))
    if tmpIBC >= .9999999:
        return True
    else:
        return False

def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())

#https://stackoverflow.com/questions/4417546/constantly-print-subprocess-output-while-process-is-running
def execute(cmd):
    popen = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line
    popen.stdout.close()
    return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)

def gen_filenames(prefix, suffix, places=3):
    """Generate sequential filenames with the format <prefix><index><suffix>

       The index field is padded with leading zeroes to the specified number of places

       http://stackoverflow.com/questions/5068461/how-do-you-increment-file-name-in-python
    """
    pattern = "{}{{:0{}d}}{}".format(prefix, places, suffix)
    for i in count(1):
        yield pattern.format(i)

def fastmech_change_cd(hdf_file, newCd):
    # hdf5_file_name = r'F:\Kootenai Project\USACE\Braided\Case11_tmp.cgn'
    # r+ adds read/write permisions to file
    file = h5py.File(hdf_file, 'r+')
    group = file['/iRIC/CalculationConditions/FM_HydAttCD/Value']
    dset = group[u' data']
    # print dset[0]
    dset[0] = newCd
    # print dset[0]
    file.close()

def fastmech_BCs(hdf_file, Q, H_DS, H_US, iniType, OneDCD):
    file = h5py.File(hdf_file, 'r+')
    group = file['/iRIC/CalculationConditions/FM_HydAttQ/Value']
    dset = group[u' data']
    dset[0] = Q
    group2 = file['/iRIC/CalculationConditions/FM_HydAttWS2/Value']
    dset2 = group2[u' data']
    dset2[0] = H_US
    group3 = file['/iRIC/CalculationConditions/FM_HydAttWS/Value']
    dset3 = group3[u' data']
    dset3[0] = H_DS
    group4 = file['/iRIC/CalculationConditions/FM_HydAttWSType/Value']
    dset4 = group4[u' data']
    dset4[0] = iniType
    # group5 = file['/iRIC/CalculationConditions/FM_HydAttWS1DStage/Value']
    # dset5 = group5[u' data']
    # dset5[0] = OneDStage
    # group6 = file['/iRIC/CalculationConditions/FM_HydAttWS1DDisch/Value']
    # dset6 = group6[u' data']
    # dset6[0] = OneDQ
    group7 = file['/iRIC/CalculationConditions/FM_HydAttWS1DCD/Value']
    dset7 = group7[u' data']
    dset7[0] = OneDCD
    file.close()

def fastmech_params(hdf_file, Itermax, endLEV):
    file = h5py.File(hdf_file, 'r+')
    group = file['/iRIC/CalculationConditions/FM_EndLEV/Value']
    dset = group[u' data']
    dset[0] = endLEV
    group = file['/iRIC/CalculationConditions/FM_SolAttItm/Value']
    dset = group[u' data']
    dset[0] = Itermax
    file.close()

def fastmech_change_var_cd(hdf_file, newCd_0, newCd_1):
    # hdf5_file_name = r'F:\Kootenai Project\USACE\Braided\Case11_tmp.cgn'
    # r+ adds read/write permisions to file
    file = h5py.File(hdf_file, 'r+')
    group = file['/iRIC/iRICZone/GridConditions/sanddepth/Value']
    dset = group[u' data']
    group2 = file['/iRIC/iRICZone/GridConditions/roughness/Value']
    dset2 = group2[u' data']
    for index, val in enumerate(dset):
        if val == 0.0:
            dset2[index] = newCd_0
        # else:
            # dset2[index] = newCd_1 #keep values in original project, change only values with 0
    # print dset[0]
    # print dset[0]
    file.close()

def add_fastmech_solver_to_path(solverPath):
    print(os.environ['PATH'])
    os.environ['PATH'] += solverPath
    print("\n")
    print('new path')
    print(os.environ['PATH'])

def add_to_path(path):
    print(os.environ['PATH'])
    os.environ['PATH'] += path
    print("\n")
    print('new path')
    print(os.environ['PATH'])

def create_vtk_structured_grid(sgrid, hdf5_file_name, xoffset, yoffset):
    # type: (object) -> object
    file = h5py.File(hdf5_file_name, 'r')
    xcoord_grp = file['/iRIC/iRICZone/GridCoordinates/CoordinateX']
    print(xcoord_grp.keys())
    ycoord_grp = file['/iRIC/iRICZone/GridCoordinates/CoordinateY']
    print(ycoord_grp.keys())
    wse_grp = file['iRIC/iRICZone/FlowSolution1/WaterSurfaceElevation']
    print(wse_grp.keys())
    topo_grp = file['iRIC/iRICZone/FlowSolution1/Elevation']
    print(topo_grp.keys())
    ibc_grp = file['iRIC/iRICZone/FlowSolution1/IBC']
    velx_grp = file['iRIC/iRICZone/FlowSolution1/VelocityX']
    vely_grp = file['iRIC/iRICZone/FlowSolution1/VelocityY']

    xcoord_data = xcoord_grp[u' data']
    ycoord_data = ycoord_grp[u' data']
    wse_data = wse_grp[u' data']
    topo_data = topo_grp[u' data']
    ibc_data = ibc_grp[u' data']
    velx_data = velx_grp[u' data']
    vely_data = vely_grp[u' data']


    # SGrid = vtk.vtkStructuredGrid()
    ny, nx, = xcoord_data.shape
    print(ny, nx)
    sgrid.SetDimensions(nx, ny, 1)
    points = vtk.vtkPoints()
    wseVal = vtk.vtkFloatArray()
    wseVal.SetNumberOfComponents(1)
    ibcVal = vtk.vtkIntArray()
    ibcVal.SetNumberOfComponents(1)
    velVal = vtk.vtkFloatArray()
    velVal.SetNumberOfComponents(1)
    for j in range(ny):
        for i in range(nx):
            points.InsertNextPoint(xcoord_data[j, i] - xoffset, ycoord_data[j, i] - yoffset, 0.0)
            wseVal.InsertNextValue(wse_data[j, i])
            ibcVal.InsertNextValue(ibc_data[j, i])
            velVal.InsertNextValue(np.sqrt(np.power(velx_data[j, i],2) + np.power(vely_data[j,i],2)))
        sgrid.SetPoints(points)

        sgrid.GetPointData().AddArray(wseVal)
        sgrid.GetPointData().AddArray(ibcVal)
        sgrid.GetPointData().AddArray(velVal)
    wseVal.SetName("WSE")
    ibcVal.SetName("IBC")
    velVal.SetName("Velocity")

print('Number of arguments:', len(sys.argv), 'arguments.')
print('Argument List:', str(sys.argv))
setfile = sys.argv[1]
config = configparser.ConfigParser()
config.read(setfile)

# add_to_path(config.get('Params','lib_path'))
add_to_path(config.get('Params', 'hdfdiff_path'))
add_fastmech_solver_to_path(config.get('Params','solver_path'))

g = gen_filenames("FM_Test_", ".cgns")

# TEST #2
hdf5_file_name1 = next(g)

copyfile(config.get('Params','base_test_1_file'),
            hdf5_file_name1)

for path in execute(["Fastmech.exe", hdf5_file_name1]):
    print('Test1',path, end="")



# TEST #2
hdf5_file_name2 = next(g)

copyfile(config.get('Params','base_test_2_file'),
            hdf5_file_name2)

for path in execute(["Fastmech.exe", hdf5_file_name2]):
    print('Test2',path, end="")




# TEST #3
hdf5_file_name3 = next(g)

copyfile(config.get('Params','base_test_3_file'),
            hdf5_file_name3)

for path in execute(["Fastmech.exe", hdf5_file_name3]):
    print('Test3',path, end="")


#SUMMARY
print('Summary')
result1 = subprocess.run(["h5diff.exe", config.get('Params','base_test_1_file'), hdf5_file_name1],
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
if result1.returncode == 0:
    print("No differences found in Test1")
else:
    result1 = subprocess.run(["h5diff.exe", config.get('Params', 'base_test_1_file'), hdf5_file_name1],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print("Differences found in Test1")
    print(result1.stdout)

result2 = subprocess.run(["h5diff.exe", config.get('Params','base_test_2_file'), hdf5_file_name2],
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
if result2.returncode == 0:
    print("No differences found in Test2")
else:
    result2 = subprocess.run(["h5diff.exe", config.get('Params', 'base_test_2_file'), hdf5_file_name2],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print("Differences found in Test2")
    print(result2.stdout)

result3 = subprocess.run(["h5diff.exe", config.get('Params','base_test_3_file'), hdf5_file_name3],
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
if result3.returncode == 0:
    print("No differences found in Test3")
else:
    result3 = subprocess.run(["h5diff.exe", config.get('Params', 'base_test_3_file'), hdf5_file_name3],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print("Differences found in Test3")
    print(result3.stdout)