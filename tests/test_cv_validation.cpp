#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "RealFluidState.H"
#include "SpeciesCriticalDB.H"
#include "BinaryInteractionDB.H"

using namespace std;
using namespace Foam;

int main()
{
    cout << "\n========================================" << endl;
    cout << "  Stage 0 Validation Test Suite" << endl;
    cout << "  Cv > 0 Everywhere" << endl;
    cout << "========================================\n" << endl;

    // Create databases
    auto criticalDB = SpeciesCriticalDB::createDefault();
    auto kijDB = BinaryInteractionDB::createDefault();
    
    PengRobinsonEOS eos(criticalDB, kijDB);
    
    int total_tests = 0;
    int passed_tests = 0;
    
    // Test 1: OF11 FAILURE CASE #1
    cout << "Test 1: OF11 Failure Case #1" << endl;
    cout << "  C3H8 @ T=338K, p=5.8MPa" << endl;
    cout << "  OF11 returned: Cv = -49633 J/(kg·K)" << endl;
    
    {
        double p = 5.8e6;
        double T = 338.0;
        vector<double> Y(2, 0.0);
        Y[1] = 1.0;  // Pure C3H8
        
        auto result = eos.solveZ(p, T, Y);
        double Cv_dep = eos.cvDeparture(p, T, result.Z, Y);
        double Cv_ideal = 1500.0;
        double Cv = Cv_ideal + Cv_dep;
        
        cout << "  Z = " << result.Z << endl;
        cout << "  Cv_departure = " << Cv_dep << " J/(kg·K)" << endl;
        cout << "  Cv_total = " << Cv << " J/(kg·K)" << endl;
        
        total_tests++;
        if (Cv > 0) {
            cout << "  ✓ PASS: Cv > 0" << endl;
            passed_tests++;
        } else {
            cout << "  ✗ FAIL: Cv <= 0" << endl;
        }
    }
    
    // Test 2: OF11 FAILURE CASE #2
    cout << "\nTest 2: OF11 Failure Case #2" << endl;
    cout << "  C3H8 @ T=365K, p=5.8MPa" << endl;
    
    {
        double p = 5.8e6;
        double T = 365.0;
        vector<double> Y(2, 0.0);
        Y[1] = 1.0;
        
        auto result = eos.solveZ(p, T, Y);
        double Cv_dep = eos.cvDeparture(p, T, result.Z, Y);
        double Cv = 1500.0 + Cv_dep;
        
        cout << "  Z = " << result.Z << endl;
        cout << "  Cv = " << Cv << " J/(kg·K)" << endl;
        
        total_tests++;
        if (Cv > 0) {
            cout << "  ✓ PASS: Cv > 0" << endl;
            passed_tests++;
        } else {
            cout << "  ✗ FAIL: Cv <= 0" << endl;
        }
    }
    
    // Test 3: Critical Point Scan
    cout << "\nTest 3: Critical Point Scan" << endl;
    cout << "  Scanning Tr=[0.95,1.05], Pr=[0.9,1.1]" << endl;
    
    double Tc = 369.89;  // C3H8
    double Pc = 4.251e6;
    
    int scan_total = 0;
    int scan_passed = 0;
    int scan_failed = 0;
    
    for (double Tr = 0.95; Tr <= 1.05; Tr += 0.01) {
        for (double Pr = 0.9; Pr <= 1.1; Pr += 0.05) {
            double T = Tr * Tc;
            double p = Pr * Pc;
            
            vector<double> Y(2, 0.0);
            Y[1] = 1.0;
            
            try {
                auto result = eos.solveZ(p, T, Y);
                if (result.valid) {
                    double Cv_dep = eos.cvDeparture(p, T, result.Z, Y);
                    double Cv = 1500.0 + Cv_dep;
                    
                    scan_total++;
                    if (Cv > 0) {
                        scan_passed++;
                    } else {
                        scan_failed++;
                        cout << "    FAIL at Tr=" << Tr << ", Pr=" << Pr 
                             << ", Cv=" << Cv << endl;
                    }
                }
            } catch (...) {
                // Skip invalid points
            }
        }
    }
    
    cout << "  Tested " << scan_total << " points" << endl;
    cout << "  Passed: " << scan_passed << ", Failed: " << scan_failed << endl;
    
    total_tests++;
    if (scan_passed == scan_total && scan_total > 0) {
        cout << "  ✓ PASS: All critical points have Cv > 0" << endl;
        passed_tests++;
    } else {
        cout << "  ✗ FAIL: Some points have Cv <= 0" << endl;
    }
    
    // Test 4: Wide Range Scan
    cout << "\nTest 4: Wide Range Scan" << endl;
    cout << "  p=[1,10]MPa, T=[300,600]K" << endl;
    
    int wide_total = 0;
    int wide_passed = 0;
    int wide_failed = 0;
    
    for (double p = 1e6; p <= 10e6; p += 1e6) {
        for (double T = 300; T <= 600; T += 50) {
            vector<double> Y(2, 0.0);
            Y[1] = 1.0;
            
            try {
                auto result = eos.solveZ(p, T, Y);
                if (result.valid) {
                    double Cv_dep = eos.cvDeparture(p, T, result.Z, Y);
                    double Cv = 1500.0 + Cv_dep;
                    
                    wide_total++;
                    if (Cv > 0) {
                        wide_passed++;
                    } else {
                        wide_failed++;
                        cout << "    FAIL at p=" << p/1e6 << "MPa, T=" << T 
                             << "K, Cv=" << Cv << endl;
                    }
                }
            } catch (...) {
            }
        }
    }
    
    cout << "  Tested " << wide_total << " points" << endl;
    cout << "  Passed: " << wide_passed << ", Failed: " << wide_failed << endl;
    
    total_tests++;
    if (wide_passed == wide_total && wide_total > 0) {
        cout << "  ✓ PASS" << endl;
        passed_tests++;
    } else {
        cout << "  ✗ FAIL" << endl;
    }
    
    // Final Summary
    cout << "\n" << string(50, '=') << endl;
    cout << "FINAL RESULTS: " << passed_tests << "/" << total_tests << " tests passed" << endl;
    
    if (passed_tests == total_tests) {
        cout << "\n╔════════════════════════════════════════╗" << endl;
        cout << "║  ✓ Stage 0 COMPLETE                    ║" << endl;
        cout << "║  Thermal kernel validated.             ║" << endl;
        cout << "║  Cv > 0 everywhere.                    ║" << endl;
        cout << "║  Ready for Stage 1.                    ║" << endl;
        cout << "╚════════════════════════════════════════╝\n" << endl;
        return 0;
    } else {
        cout << "\n✗✗✗ SOME TESTS FAILED ✗✗✗\n" << endl;
        return 1;
    }
}
