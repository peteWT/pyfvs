"""
Native FVS Comparison Example

Demonstrates running growth simulations through both pyfvs (Python)
and the native FVS Fortran library, then comparing results.

Prerequisites:
    - Build and install the FVS shared library (see native/BUILD.md)
    - Set FVS_LIB_PATH if not installed to a standard location

Usage:
    uv run python examples/native_fvs_comparison.py
"""

from rich.console import Console
from rich.table import Table
from rich.panel import Panel

from pyfvs import Stand
from pyfvs.native import fvs_library_available, get_library_info

console = Console()


def check_library_availability():
    """Check which FVS variant libraries are available."""
    console.print(Panel("[bold]FVS Library Availability[/bold]"))

    variants = ["SN", "LS", "PN", "WC", "NE", "CS", "OP", "CA", "OC", "WS"]

    table = Table()
    table.add_column("Variant", style="bold")
    table.add_column("Available")
    table.add_column("Path")

    any_available = False
    for variant in variants:
        info = get_library_info(variant)
        available = info["available"]
        path = info["path"] or "-"
        status = "[green]Yes[/green]" if available else "[red]No[/red]"
        table.add_row(variant, status, path)
        if available:
            any_available = True

    console.print(table)
    console.print()

    if not any_available:
        console.print(
            "[yellow]No FVS shared libraries found.[/yellow]\n"
            "To build and install, see: src/pyfvs/native/BUILD.md\n"
            "Or set: export FVS_LIB_PATH=/path/to/libs\n"
        )

    return any_available


def run_pyfvs_only_demo():
    """Run a pyfvs-only simulation to demonstrate the API."""
    console.print(Panel("[bold]PyFVS Simulation (Python-only)[/bold]"))

    stand = Stand.initialize_planted(
        trees_per_acre=500,
        site_index=70,
        species="LP",
        variant="SN",
        ecounit="M231",
    )

    table = Table(title="Loblolly Pine Yield Table (SN variant, SI=70, M231)")
    table.add_column("Year", justify="center")
    table.add_column("TPA", justify="right")
    table.add_column("BA", justify="right")
    table.add_column("QMD", justify="right")
    table.add_column("Volume", justify="right")

    for year in range(5, 55, 5):
        stand.grow(years=5)
        m = stand.get_metrics()
        table.add_row(
            str(year),
            f"{m['tpa']:.0f}",
            f"{m['basal_area']:.1f}",
            f"{m['qmd']:.2f}",
            f"{m['volume']:.0f}",
        )

    console.print(table)
    console.print()


def run_native_comparison():
    """Run side-by-side comparison if a native library is available."""
    # Find the first available variant
    for variant, species, si, tpa in [
        ("SN", "LP", 70, 500),
        ("PN", "DF", 120, 400),
        ("LS", "RN", 65, 500),
    ]:
        if fvs_library_available(variant):
            console.print(
                Panel(f"[bold]PyFVS vs Native FVS: {variant} variant[/bold]")
            )

            from validation.native.compare_native import (
                compare_stand_growth,
                generate_validation_report,
            )

            results = compare_stand_growth(
                trees_per_acre=tpa,
                site_index=si,
                species=species,
                variant=variant,
                years=50,
            )
            generate_validation_report(results)
            return True

    return False


def main():
    console.print()
    console.rule("[bold blue]PyFVS Native FVS Comparison[/bold blue]")
    console.print()

    # Step 1: Check library availability
    libraries_found = check_library_availability()

    # Step 2: Always show pyfvs-only demo
    run_pyfvs_only_demo()

    # Step 3: Run comparison if libraries available
    if libraries_found:
        ran_comparison = run_native_comparison()
        if ran_comparison:
            console.print(
                "[green]Comparison complete. "
                "See validation/results/native/ for detailed results.[/green]"
            )
    else:
        console.print(
            "[dim]Install FVS shared libraries to enable native comparison.[/dim]"
        )


if __name__ == "__main__":
    main()
