"""Run a small, self-contained displacement-speed sweep."""

from wigley import make_wigley_hull
import michell


def main():
    hull = make_wigley_hull(michell)
    solver = michell.Resistance(hull)
    solver.set_grid_precision(michell.ResistanceGridPrecision.Fine)
    solver.prepare()

    print(f"GenericMichell {solver.get_version()}")
    print(f"Lwl={hull.get_length_wl():.4f} m  "
          f"Bwl={hull.get_beam_wl():.4f} m  T={hull.get_draft():.4f} m")
    print(" Fn_L    speed_m/s    R_total_N   R_friction_N   R_wave_N   R_pressure_N")
    for fn in (0.20, 0.25, 0.30, 0.35, 0.40):
        speed = hull.get_speed(fn)
        total = solver.get_total_resistance(speed, michell.ResistanceMethod.Auto)
        friction = solver.get_friction_resistance(speed, michell.ResistanceMethod.Auto)
        wave = solver.get_wave_resistance(speed, michell.ResistanceMethod.Auto)
        pressure = solver.get_pressure_resistance(speed, michell.ResistanceMethod.Auto)
        print(f"{fn:5.2f}  {speed:11.5f}  {total:11.5f}  {friction:13.5f}  "
              f"{wave:9.5f}  {pressure:13.5f}")


if __name__ == "__main__":
    main()
