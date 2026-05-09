"""
extract_validation.py  (rev 3 - node-based, version-agnostic)
-------------------------------------------------------------
Run with ABAQUS Python:

    abaqus python extract_validation.py <jobname>

Reads only NODES and ELEMENT CONNECTIVITY from the ODB. These are written by
every ABAQUS job regardless of step type, version, or field-output settings,
and so this script is robust against the OdbInstance.getMassProperties() and
field-output naming differences seen across ABAQUS releases.

For each CPS4R element (axis-aligned quadrilateral on our Cartesian mesh),
we compute:
    area_e   = dx * dy  via the shoelace formula
    centre_e = mean of 4 node coordinates  (exact for parallelograms)

Then accumulate to area, centroid, and centroidal area-moments of inertia.

Verbose stdout reports each step so any failure point is visible in the
batch log.
"""

from __future__ import print_function

import os
import sys
import traceback

from odbAccess import openOdb


def shoelace_area(xs, ys):
    n = len(xs)
    s = 0.0
    for i in range(n):
        j = (i + 1) % n
        s += xs[i] * ys[j] - xs[j] * ys[i]
    return abs(s) * 0.5


def main():
    if len(sys.argv) < 2:
        print("usage: abaqus python extract_validation.py <jobname>")
        sys.exit(1)
    name = sys.argv[1]
    odb_path = name if name.endswith(".odb") else name + ".odb"
    job_name = odb_path.replace(".odb", "")

    print("=" * 60)
    print("Job: " + job_name)
    print("=" * 60)

    if not os.path.exists(odb_path):
        print("ERROR: missing odb: " + odb_path)
        sys.exit(2)

    odb = openOdb(odb_path, readOnly=True)
    try:
        # 1. Inspect what is in the ODB
        instance_names = list(odb.rootAssembly.instances.keys())
        print("Instances in odb: " + str(instance_names))
        if not instance_names:
            print("ERROR: no instance in ODB")
            sys.exit(3)
        inst_name = instance_names[0]
        inst = odb.rootAssembly.instances[inst_name]
        n_nodes = len(inst.nodes)
        n_elems = len(inst.elements)
        print("Instance '" + inst_name + "' nodes=" + str(n_nodes) +
              "  elements=" + str(n_elems))

        # 2. Build node lookup table:  label -> (x, y)
        nodes = {}
        for nd in inst.nodes:
            c = nd.coordinates
            nodes[nd.label] = (float(c[0]), float(c[1]))
        print("Loaded " + str(len(nodes)) + " nodes")

        # 3. First pass: total area and first moments
        A = 0.0
        sx = 0.0
        sy = 0.0
        elem_data = []                 # cache for the second pass
        for el in inst.elements:
            conn = el.connectivity     # tuple of 4 node labels for CPS4R
            xs = [nodes[i][0] for i in conn]
            ys = [nodes[i][1] for i in conn]
            a = shoelace_area(xs, ys)
            xc = sum(xs) / float(len(xs))
            yc = sum(ys) / float(len(ys))
            A += a
            sx += xc * a
            sy += yc * a
            elem_data.append((a, xc, yc))
        x_c = sx / A
        y_c = sy / A
        print("Pass 1:  A=" + str(round(A, 6)) +
              "  centroid=(" + str(round(x_c, 6)) + ", "
              + str(round(y_c, 6)) + ")")

        # 4. Second pass: second moments about centroid
        Ixx = 0.0
        Iyy = 0.0
        Ixy = 0.0
        for a, xc, yc in elem_data:
            dx = xc - x_c
            dy = yc - y_c
            Ixx += dy * dy * a
            Iyy += dx * dx * a
            Ixy += dx * dy * a
        print("Pass 2:  Ixx=" + str(round(Ixx, 6)) +
              "  Iyy=" + str(round(Iyy, 6)) +
              "  Ixy=" + str(round(Ixy, 6)))

        # 5. Append to CSV
        out_csv = os.path.join(os.path.dirname(odb_path) or ".",
                               "validation_results.csv")
        write_header = not os.path.exists(out_csv)
        f = open(out_csv, "a")
        try:
            if write_header:
                f.write("job,n_elems,area,x_c,y_c,I_xx,I_yy,I_xy\n")
            f.write("{job},{n},{A:.10f},{x:.10f},{y:.10f},"
                    "{Ixx:.10f},{Iyy:.10f},{Ixy:.10f}\n".format(
                        job=job_name, n=n_elems,
                        A=A, x=x_c, y=y_c,
                        Ixx=Ixx, Iyy=Iyy, Ixy=Ixy))
        finally:
            f.close()
        print("WROTE row to " + out_csv)
    except Exception:
        print("UNHANDLED EXCEPTION while processing " + job_name + ":")
        traceback.print_exc()
        sys.exit(5)
    finally:
        odb.close()


if __name__ == "__main__":
    main()
