#!/usr/bin/env python3
"""Create the specimen and heater CAD from config/geometry.yaml."""

import math
import sys
from pathlib import Path

import FreeCAD as App
import Mesh
import Part
import yaml


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config" / "geometry.yaml"
OUT = ROOT / "cad" / "generated"


def point(radius: float, angle_rad: float, z: float = 0.0) -> App.Vector:
    return App.Vector(radius * math.cos(angle_rad), radius * math.sin(angle_rad), z)


def heater_sector(inner_r: float, outer_r: float, z0: float, height: float,
                  center_angle: float, span: float) -> Part.Shape:
    a0 = center_angle - span / 2.0
    am = center_angle
    a1 = center_angle + span / 2.0

    outer = Part.Arc(point(outer_r, a0), point(outer_r, am), point(outer_r, a1)).toShape()
    side1 = Part.makeLine(point(outer_r, a1), point(inner_r, a1))
    inner = Part.Arc(point(inner_r, a1), point(inner_r, am), point(inner_r, a0)).toShape()
    side0 = Part.makeLine(point(inner_r, a0), point(outer_r, a0))
    face = Part.Face(Part.Wire([outer, side1, inner, side0]))
    shape = face.extrude(App.Vector(0, 0, height))
    shape.translate(App.Vector(0, 0, z0))
    return shape


def main() -> None:
    cfg = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    specimen_cfg = cfg["specimen"]
    heater_cfg = cfg["heater"]

    outer_r = specimen_cfg["outer_diameter_mm"] / 2.0
    inner_r = specimen_cfg["inner_diameter_mm"] / 2.0
    height = specimen_cfg["height_mm"]
    heater_height = heater_cfg["axial_height_mm"]
    heater_z0 = heater_cfg["center_height_mm"] - heater_height / 2.0
    heater_outer_r = outer_r + heater_cfg["thickness_mm"]
    heater_span = heater_cfg["circumferential_length_mm"] / outer_r
    heater_center = math.radians(heater_cfg["center_angle_deg"])

    specimen = Part.makeCylinder(outer_r, height).cut(Part.makeCylinder(inner_r, height))
    heater = heater_sector(
        outer_r, heater_outer_r, heater_z0, heater_height, heater_center, heater_span
    )
    # Numerical selection volume. The 0.5 mm overlap/extension is a topoSet
    # tolerance, not a physical heater dimension.
    heater_selection = heater_sector(
        outer_r - 0.5,
        heater_outer_r + 0.5,
        heater_z0 - 0.5,
        heater_height + 1.0,
        heater_center,
        heater_span,
    )

    OUT.mkdir(parents=True, exist_ok=True)
    doc = App.newDocument("hollow_cylinder_heater")
    specimen_obj = doc.addObject("PartDesign::Feature", "Specimen")
    specimen_obj.Label = "Carbon steel hollow cylinder"
    specimen_obj.Shape = specimen
    heater_obj = doc.addObject("PartDesign::Feature", "HeaterMat")
    heater_obj.Label = "Silicone heater mat"
    heater_obj.Shape = heater

    params = doc.addObject("App::FeaturePython", "Parameters")
    for name, value in {
        "OuterDiameter": specimen_cfg["outer_diameter_mm"],
        "InnerDiameter": specimen_cfg["inner_diameter_mm"],
        "SpecimenHeight": height,
        "HeaterHeight": heater_height,
        "HeaterArcLength": heater_cfg["circumferential_length_mm"],
        "HeaterThickness": heater_cfg["thickness_mm"],
    }.items():
        params.addProperty("App::PropertyLength", name)
        setattr(params, name, value)

    doc.recompute()
    doc.saveAs(str(OUT / "hollow_cylinder_heater.FCStd"))
    Part.export([specimen_obj, heater_obj], str(OUT / "hollow_cylinder_heater.step"))
    Mesh.export([specimen_obj], str(OUT / "specimen_mm.stl"))
    Mesh.export([heater_obj], str(OUT / "heaterMat_mm.stl"))
    selection_obj = doc.addObject("PartDesign::Feature", "HeaterSelection")
    selection_obj.Label = "Numerical face selection volume"
    selection_obj.Shape = heater_selection
    doc.recompute()
    doc.saveAs(str(OUT / "hollow_cylinder_heater.FCStd"))
    Mesh.export([selection_obj], str(OUT / "heaterSelection_mm.stl"))

    volume_cm3 = specimen.Volume / 1000.0
    density = specimen_cfg["measured_mass_kg"] / (specimen.Volume * 1.0e-9)
    print(f"specimen volume: {volume_cm3:.3f} cm3")
    print(f"density from measured mass: {density:.1f} kg/m3")
    print(f"heater angular span: {math.degrees(heater_span):.3f} deg")
    print(f"outputs: {OUT}")


# FreeCADCmd executes a script as a module, so __name__ is not guaranteed to
# be "__main__". The file is an executable macro and intentionally runs here.
main()
