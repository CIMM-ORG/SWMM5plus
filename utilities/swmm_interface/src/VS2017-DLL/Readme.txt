INSTRUCTIONS FOR COMPILING SWMM5.DLL USING MICROSOFT VISUAL STUDIO 2017
=======================================================================

1. Create a sub-directory named VS2017-DLL under the directory where
   the SWMM 5 Engine source code files are stored and copy SWMM5.DEF
   and VS2017-DLL.VCPROJ to it.

2. Launch Visual Studio 2017 and use the File >> Open command to open
   the VS2017-DLL.VCPROJ file.

3. Issue the Build >> Configuration Manager command and select the
   Release configuration.

4. Issue the Build >> Build VS2017-DLL command to build SWMM5.DLL
   (which will appear in the Release subdirectory underneath the
   VC2017-DLL directory).
