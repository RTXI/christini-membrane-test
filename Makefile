PLUGIN_NAME = membrane_test

HEADERS = Membrane_Test.h \
	Membrane_Test_MainWindow_UI.h

SOURCES = Membrane_Test.cpp moc_Membrane_Test.cpp

LIBS = -lgsl -lgslcblas -lrtmath

### Do not edit below this line ###

include $(shell rtxi_plugin_config --pkgdata-dir)/Makefile.plugin_compile

Membrane_Test.cpp: Membrane_Test_MainWindow_UI.h

Membrane_Test_MainWindow_UI.h: Membrane_Test_MainWindow.ui
	uic Membrane_Test_MainWindow.ui -o Membrane_Test_MainWindow_UI.h

clean: uiclean

uiclean:
	rm Membrane_Test_MainWindow_UI.h
