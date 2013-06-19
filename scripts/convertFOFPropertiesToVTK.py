import vtk
import os
import os.path
import subprocess
import re
import sys
import StringIO

def main(inputfile, outputfile):
    GENERICIOBIN = r'/home/casey.goodlett/projects/sciviz/generic-bin'
    INSPECTOR = os.path.join(GENERICIOBIN, 'GenericIOInspector')

    output = subprocess.check_output([INSPECTOR, '--file', inputfile, '--all'])

    elementsmatcher = re.compile(r'NUMBER OF ELEMENTS=([0-9]+)')

    points = vtk.vtkPoints()

    count = vtk.vtkIntArray()
    count.SetNumberOfComponents(1)
    count.SetName('count')

    tag = vtk.vtkIdTypeArray()
    tag.SetNumberOfComponents(1)
    tag.SetName('tag')

    mass = vtk.vtkDoubleArray()
    mass.SetNumberOfComponents(1)
    mass.SetName('mass')

    meanPosition = vtk.vtkDoubleArray()
    meanPosition.SetNumberOfComponents(3)
    meanPosition.SetName('mean')

    vel = vtk.vtkDoubleArray()
    vel.SetNumberOfComponents(3)
    vel.SetName('velocity')

    disp = vtk.vtkDoubleArray()
    disp.SetNumberOfComponents(1)
    disp.SetName('velocity_dispersal')

    arrays = [count, tag, mass, meanPosition, vel, disp]

    i = 0
    for line in StringIO.StringIO(output):
        if line[0] == '#':
            result = elementsmatcher.search(line)
            if result:
                numberElements = int(result.group(1))
                points.SetNumberOfPoints(numberElements)
                for array in arrays:
                    array.SetNumberOfTuples(numberElements)
        else:
            data = line.split()
            c, t, m, cx, cy, cz, mx, my, mz, vx, vy, vz, veldisp = data

            count.SetTuple1(i, int(c))
            tag.SetTuple1(i, int(t))
            mass.SetTuple1(i, float(m))
            points.SetPoint(i, [float(x) for x in [cx, cy, cz]])
            meanPosition.SetTuple(i, [float(x) for x in [mx, my, mz]])
            vel.SetTuple(i, [float(x) for x in [vx, vy, vz]])
            disp.SetTuple1(i, float(veldisp))

            i = i + 1

    pd = vtk.vtkPolyData()
    pd.SetPoints(points)
    for array in arrays:
        pd.GetPointData().AddArray(array)

    mask = vtk.vtkMaskPoints()
    mask.SetOnRatio(1)
    mask.RandomModeOff()
    mask.GenerateVerticesOn()
    mask.SingleVertexPerCellOn()

    # Change this for VTK 6
    mask.SetInput(pd)

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetInputConnection(mask.GetOutputPort())
    writer.SetFileName(outputfile)
    writer.Update()

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
