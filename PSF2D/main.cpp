#include "mainwindow.h"
#include <QApplication>
#include <viennacl/vector.hpp>
#include <viennacl/matrix.hpp>
#include <viennacl/ocl/device.hpp>
#include <viennacl/ocl/platform.hpp>
#include <viennacl/device_specific/builtin_database/common.hpp>
#include <Eigen/Eigen>
#include <iostream>

int main(int argc, char *argv[])
{
    typedef std::vector< viennacl::ocl::platform > platforms_type;
    platforms_type platforms = viennacl::ocl::get_platforms();

    bool is_first_element = true;
    for (platforms_type::iterator platform_iter  = platforms.begin();
                                  platform_iter != platforms.end();
                                  ++platform_iter)
    {
        typedef std::vector<viennacl::ocl::device> devices_type;
        devices_type devices = platform_iter->devices(CL_DEVICE_TYPE_ALL);

        std::cout << "# =========================================" << std::endl;
        std::cout << "#         Platform Information             " << std::endl;
        std::cout << "# =========================================" << std::endl;

        std::cout << "#" << std::endl;
        std::cout << "# Vendor and version: " << platform_iter->info() << std::endl;
        std::cout << "#" << std::endl;

        if (is_first_element)
        {
            std::cout << "# ViennaCL uses this OpenCL platform by default." << std::endl;
            is_first_element = false;
        }

        std::cout << "# " << std::endl;
        std::cout << "# Available Devices: " << std::endl;
        std::cout << "# " << std::endl;
        for (devices_type::iterator iter = devices.begin(); iter != devices.end(); iter++)
        {
            std::cout << std::endl;

            std::cout << "  -----------------------------------------" << std::endl;
            std::cout << iter->full_info();
            std::cout << "ViennaCL Device Architecture:  " << iter->architecture_family() << std::endl;
            std::cout << "ViennaCL Database Mapped Name: " << viennacl::device_specific::builtin_database::get_mapped_device_name(iter->name(), iter->vendor_id()) << std::endl;
            std::cout << "  -----------------------------------------" << std::endl;
        }
        std::cout << std::endl;
        std::cout << "###########################################" << std::endl;
        std::cout << std::endl;
    }
    Eigen::MatrixXd eigenMat = Eigen::MatrixXd::Random(10,10);
    Eigen::VectorXd eigenVec = Eigen::VectorXd::Random(10);
    viennacl::vector<double> vclVector(10);
    viennacl::vector<double> vclResult(10);
    viennacl::matrix<double> vclDenseMatrix(10,10);
    viennacl::copy(eigenVec,vclVector);
    viennacl::copy(eigenMat,vclDenseMatrix);
    vclResult = viennacl::linalg::prod(vclDenseMatrix,vclVector);
    std::cout<<vclResult<<std::endl;
    QApplication a(argc, argv);
    MainWindow w;
    w.show();

    return a.exec();
}
