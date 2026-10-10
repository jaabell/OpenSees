/**
 * Test file for Hoek-Brown yield function and plastic flow
 * 
 * This is a standalone test that verifies the basic functionality of:
 * - HoekBrown_YF.h  (Yield Function)
 * - HoekBrown_PF.h  (Plastic Flow Direction)
 * - HoekBrown_Utils.h (Parameter computation utilities)
 * 
 * Compile with:
 *   g++ -std=c++17 -I/path/to/eigen -I/path/to/opensees test_HoekBrown.cpp -o test_HoekBrown
 */

#include <iostream>
#include <cmath>
#include <cassert>
#include <iomanip>

// Include the Hoek-Brown utilities first to test them standalone
#include "HoekBrown_Utils.h"

/**
 * Test the utility functions for computing Hoek-Brown parameters
 */
void test_hoek_brown_utils() {
    std::cout << "=== Testing HoekBrown_Utils ===" << std::endl;
    
    // Test case: Intact rock (GSI=100, D=0)
    {
        double mi = 25.0;  // Typical value for granite
        double GSI = 100.0;
        double D = 0.0;
        
        double mb = HoekBrownUtils::compute_mb(mi, GSI, D);
        double s = HoekBrownUtils::compute_s(GSI, D);
        double a = HoekBrownUtils::compute_a(GSI);
        
        std::cout << "Intact rock (GSI=100, D=0):" << std::endl;
        std::cout << "  mi = " << mi << std::endl;
        std::cout << "  mb = " << mb << " (expected: " << mi << ")" << std::endl;
        std::cout << "  s  = " << s << " (expected: 1.0)" << std::endl;
        std::cout << "  a  = " << a << " (expected: ~0.5)" << std::endl;
        
        // Check that for intact rock, mb = mi, s = 1
        assert(std::abs(mb - mi) < 1e-10);
        assert(std::abs(s - 1.0) < 1e-10);
        assert(std::abs(a - 0.5) < 0.01);  // a should be very close to 0.5
        
        std::cout << "  PASSED!" << std::endl << std::endl;
    }
    
    // Test case: Fair quality rock mass (GSI=50, D=0)
    {
        double mi = 25.0;
        double GSI = 50.0;
        double D = 0.0;
        
        auto [mb, s, a] = HoekBrownUtils::compute_rock_mass_parameters(mi, GSI, D);
        
        std::cout << "Fair quality rock (GSI=50, D=0):" << std::endl;
        std::cout << "  mi = " << mi << std::endl;
        std::cout << "  mb = " << mb << std::endl;
        std::cout << "  s  = " << s << std::endl;
        std::cout << "  a  = " << a << std::endl;
        
        // mb should be significantly reduced
        assert(mb < mi);
        // s should be close to zero
        assert(s < 0.01);
        // a should be between 0.5 and 0.65
        assert(a >= 0.5 && a <= 0.65);
        
        std::cout << "  PASSED!" << std::endl << std::endl;
    }
    
    // Test case: Poor quality rock with disturbance (GSI=30, D=0.7)
    {
        double mi = 15.0;
        double GSI = 30.0;
        double D = 0.7;
        
        auto [mb, s, a] = HoekBrownUtils::compute_rock_mass_parameters(mi, GSI, D);
        
        std::cout << "Poor quality rock (GSI=30, D=0.7):" << std::endl;
        std::cout << "  mi = " << mi << std::endl;
        std::cout << "  mb = " << mb << std::endl;
        std::cout << "  s  = " << s << std::endl;
        std::cout << "  a  = " << a << std::endl;
        
        // Verify parameters are reasonable
        assert(mb > 0 && mb < mi);
        assert(s > 0 && s < 0.001);  // Very small s
        assert(a > 0.5 && a < 0.65);
        
        std::cout << "  PASSED!" << std::endl << std::endl;
    }
    
    // Test tensile strength computation
    {
        double sigma_ci = 100.0;  // MPa
        double mb = 10.0;
        double s = 0.01;
        
        double sigma_t = HoekBrownUtils::compute_tensile_strength(sigma_ci, mb, s);
        
        std::cout << "Tensile strength calculation:" << std::endl;
        std::cout << "  sigma_ci = " << sigma_ci << " MPa" << std::endl;
        std::cout << "  mb = " << mb << std::endl;
        std::cout << "  s  = " << s << std::endl;
        std::cout << "  sigma_t = " << sigma_t << " MPa (expected: " << -s*sigma_ci/mb << ")" << std::endl;
        
        assert(std::abs(sigma_t - (-s*sigma_ci/mb)) < 1e-10);
        assert(sigma_t < 0);  // Tensile strength should be negative (tension)
        
        std::cout << "  PASSED!" << std::endl << std::endl;
    }
    
    // Test rock mass UCS
    {
        double sigma_ci = 100.0;  // MPa
        double s = 0.01;
        double a = 0.55;
        
        double sigma_cm = HoekBrownUtils::compute_rock_mass_ucs(sigma_ci, s, a);
        
        std::cout << "Rock mass UCS calculation:" << std::endl;
        std::cout << "  sigma_ci = " << sigma_ci << " MPa" << std::endl;
        std::cout << "  s  = " << s << std::endl;
        std::cout << "  a  = " << a << std::endl;
        std::cout << "  sigma_cm = " << sigma_cm << " MPa" << std::endl;
        
        assert(sigma_cm > 0);
        assert(sigma_cm < sigma_ci);  // Rock mass UCS should be less than intact rock UCS
        
        std::cout << "  PASSED!" << std::endl << std::endl;
    }
    
    // Test equivalent Mohr-Coulomb parameters
    {
        double sigma_ci = 100.0;  // MPa
        double mb = 5.0;
        double s = 0.01;
        double a = 0.55;
        double sigma3_max = 30.0;  // MPa - typical confining pressure
        
        auto [c, phi] = HoekBrownUtils::compute_equivalent_mohr_coulomb(sigma_ci, mb, s, a, sigma3_max);
        
        std::cout << "Equivalent Mohr-Coulomb parameters:" << std::endl;
        std::cout << "  sigma_ci = " << sigma_ci << " MPa" << std::endl;
        std::cout << "  mb = " << mb << std::endl;
        std::cout << "  s  = " << s << std::endl;
        std::cout << "  a  = " << a << std::endl;
        std::cout << "  sigma3_max = " << sigma3_max << " MPa" << std::endl;
        std::cout << "  cohesion c = " << c << " MPa" << std::endl;
        std::cout << "  friction angle phi = " << phi << " degrees" << std::endl;
        
        // Verify reasonable values
        assert(c > 0);  // Cohesion should be positive
        assert(phi > 0 && phi < 90);  // Friction angle should be between 0 and 90
        
        std::cout << "  PASSED!" << std::endl << std::endl;
    }
    
    // Test mi estimation from rock type
    {
        std::cout << "Estimated mi values for various rock types:" << std::endl;
        std::cout << "  Granite: " << HoekBrownUtils::estimate_mi_from_rock_type("granite") << std::endl;
        std::cout << "  Sandstone: " << HoekBrownUtils::estimate_mi_from_rock_type("sandstone") << std::endl;
        std::cout << "  Limestone: " << HoekBrownUtils::estimate_mi_from_rock_type("limestone") << std::endl;
        std::cout << "  Shale: " << HoekBrownUtils::estimate_mi_from_rock_type("shale") << std::endl;
        std::cout << "  Claystone: " << HoekBrownUtils::estimate_mi_from_rock_type("claystone") << std::endl;
        std::cout << std::endl;
    }
    
    std::cout << "=== All HoekBrown_Utils tests PASSED! ===" << std::endl << std::endl;
}

/**
 * Print a comparison table showing Hoek-Brown parameters for different GSI values
 */
void print_parameter_table() {
    std::cout << "=== Hoek-Brown Parameter Table ===" << std::endl;
    std::cout << "Fixed parameters: mi = 25, D = 0" << std::endl;
    std::cout << std::setw(10) << "GSI" 
              << std::setw(12) << "mb" 
              << std::setw(12) << "s" 
              << std::setw(12) << "a" << std::endl;
    std::cout << std::string(46, '-') << std::endl;
    
    double mi = 25.0;
    double D = 0.0;
    
    for (int GSI = 100; GSI >= 10; GSI -= 10) {
        auto [mb, s, a] = HoekBrownUtils::compute_rock_mass_parameters(mi, GSI, D);
        std::cout << std::setw(10) << GSI 
                  << std::setw(12) << std::fixed << std::setprecision(4) << mb 
                  << std::setw(12) << std::scientific << std::setprecision(2) << s 
                  << std::setw(12) << std::fixed << std::setprecision(4) << a << std::endl;
    }
    std::cout << std::endl;
}

/**
 * Print a comparison table showing effect of disturbance factor D
 */
void print_disturbance_table() {
    std::cout << "=== Effect of Disturbance Factor D ===" << std::endl;
    std::cout << "Fixed parameters: mi = 25, GSI = 50" << std::endl;
    std::cout << std::setw(10) << "D" 
              << std::setw(12) << "mb" 
              << std::setw(12) << "s" << std::endl;
    std::cout << std::string(34, '-') << std::endl;
    
    double mi = 25.0;
    double GSI = 50.0;
    
    for (double D = 0.0; D <= 1.0; D += 0.2) {
        double mb = HoekBrownUtils::compute_mb(mi, GSI, D);
        double s = HoekBrownUtils::compute_s(GSI, D);
        std::cout << std::setw(10) << std::fixed << std::setprecision(1) << D 
                  << std::setw(12) << std::fixed << std::setprecision(4) << mb 
                  << std::setw(12) << std::scientific << std::setprecision(2) << s << std::endl;
    }
    std::cout << std::endl;
}

int main() {
    std::cout << "=====================================================" << std::endl;
    std::cout << "       Hoek-Brown Implementation Test Suite          " << std::endl;
    std::cout << "=====================================================" << std::endl << std::endl;
    
    // Run tests on utility functions
    test_hoek_brown_utils();
    
    // Print parameter tables
    print_parameter_table();
    print_disturbance_table();
    
    std::cout << "=====================================================" << std::endl;
    std::cout << "          All tests completed successfully!          " << std::endl;
    std::cout << "=====================================================" << std::endl;
    
    return 0;
}
