!include "MUI.nsh"

  ;Name and file
  Name "Minidjvu 0.9m07"
  OutFile "minidjvu-mod-0.9m07.exe"

  ;Default installation folder
  InstallDir "$PROGRAMFILES\minidjvu-mod"
  
  ;Get installation folder from registry if available
  InstallDirRegKey HKCU "Software\minidjvu-mod" ""


;--------------------------------
;Interface Settings

  !define MUI_ABORTWARNING

;--------------------------------
;Pages

  !insertmacro MUI_PAGE_WELCOME
  !insertmacro MUI_PAGE_LICENSE "COPYING"
  !insertmacro MUI_PAGE_DIRECTORY
  !insertmacro MUI_PAGE_INSTFILES
  !insertmacro MUI_PAGE_FINISH
  
  !insertmacro MUI_UNPAGE_CONFIRM
  !insertmacro MUI_UNPAGE_INSTFILES
  
;--------------------------------
;Languages
 
  !insertmacro MUI_LANGUAGE "English"
  !insertmacro MUI_LANGUAGE "Russian"

;--------------------------------
;Installer Sections

Section "Install"

  SetOutPath "$INSTDIR"
  
  Delete "$INSTDIR\bash.dll"
  Delete "$INSTDIR\bash.exe"
  Delete "$INSTDIR\djvm.exe"
  Delete "$INSTDIR\jpeg_lic.html"
  Delete "$INSTDIR\jpeg62.dll"
  Delete "$INSTDIR\libtiff3.dll"
  Delete "$INSTDIR\rm.exe"
  Delete "$INSTDIR\tiff_lic.html"
  Delete "$INSTDIR\tiffinfo.exe"
  Delete "$INSTDIR\tiffsplit.exe"
  Delete "$INSTDIR\windjview.exe"
  Delete "$INSTDIR\zlib1.dll"
  Delete "$INSTDIR\README.html"
  Delete "$SMPROGRAMS\minidjvu-mod\README.lnk"
  Delete "$INSTDIR\minidjvu-mod.nsi"
  
  File "COPYING"
  File "News"
  File "bin\minidjvu-mod.exe"
  CreateDirectory $INSTDIR\doc

  CreateDirectory $SMPROGRAMS\minidjvu-mod
  CreateShortcut  $SMPROGRAMS\minidjvu-mod\Uninstall.lnk $INSTDIR\Uninstall.exe

  ;Store installation folder
  WriteRegStr HKCU "Software\minidjvu-mod" "" $INSTDIR
  
  ;Create uninstaller
  WriteUninstaller "$INSTDIR\Uninstall.exe"

SectionEnd

;Uninstaller Section

Section "Uninstall"

  Delete "$INSTDIR\COPYING"
  Delete "$INSTDIR\minidjvu-mod.exe"

  Delete "$INSTDIR\Uninstall.exe"

  RMDir "$INSTDIR\doc"
  RMDir "$INSTDIR"

  Delete "$SMPROGRAMS\minidjvu-mod\Uninstall.lnk"
  RMDir  "$SMPROGRAMS\minidjvu-mod"

  DeleteRegKey /ifempty HKCU "Software\minidjvu-mod"
SectionEnd
