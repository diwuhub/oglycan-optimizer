"""
Example: run the optimizer on the Etanercept site catalog.

This script demonstrates the complete workflow:
1. Load the O-glycosite catalog from JSON
2. Load acquisition parameters from config
3. Evaluate glycoform feasibility across all sites
"""

from oglycan.core import evaluate, load_params, load_site_catalog


def main():
    # Load the Etanercept O-glycosite catalog
    catalog = load_site_catalog("sites/etanercept.json")

    # Load default acquisition parameters (MS instrumentation settings)
    params = load_params("examples/default_acquisition_params.json")

    # Run the evaluation
    result = evaluate(catalog, params)

    # Display results
    print(f"Composite feasibility score: {result['composite_score']:.3f}")
    print(f"\nPer-model scores:")
    for model_name, score in result["sub_model_scores"].items():
        print(f"  {model_name}: {score:.3f}")

    localization_confidence = result["composite_breakdown"]["localization_confidence"]
    print(f"\nLocalization confidence: {localization_confidence:.2%}")
    print(f"Total O-glycosites evaluated: {result['total_sites']}")


if __name__ == "__main__":
    main()
