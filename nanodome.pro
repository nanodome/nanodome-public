TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS += -std=c++11 -O2 -fopenmp
QMAKE_LFLAGS += -fopenmp
LIBS += -lstdc++fs

SOURCES += main.cpp \
    species.cpp \
    particle.cpp \
    utilities.cpp \
    virtualbond.cpp \
    cnt.cpp \
    ndm_random.cpp \
    clock.cpp \
    particlebond.cpp \
    binarykernel.cpp \
    fmbinarykernel.cpp \
    majorantkernel.cpp \
    fmmajorantkernel.cpp \
    tetra_bond.cpp \
    tetra_constrainer.cpp \
    tetra_constraint.cpp \
    tetra_face.cpp \
    tetra_link.cpp \
    tetra_tetrahedron.cpp \
    tetra_vertex.cpp \
    tinyxml2/tinyxml2.cpp \
    i_o.cpp \
    dynamicpoint.cpp \
    JSONfile.cpp \
    streamline.cpp \
    streamlinefileJSON.cpp \
    streamlinefileXML.cpp \
    JSON/JSON.cpp \
    JSON/JSONValue.cpp \
    XMLfile.cpp \
    aggregate/aggregate.cpp \
    aggregate/ellipsoidaggregate.cpp \
    aggregate/fullrigidbodyaggregate.cpp \
    aggregate/pbmaggregate.cpp \
    aggregate/rattleaggregate.cpp \
    aggregate/rigidbodyaggregate.cpp \
    aggregate/spatialaggregate.cpp \
    gasphase/gasphase.cpp \
    gasphase/gasphasecv.cpp \
    gasphase/gasphaselink.cpp \
    moments/momentmodelfriedlander.cpp \
    moments/momentmodelpratsinis.cpp \
    particlephase/cgmdparticlephase.cpp \
    particlephase/constrainedlangevinparticlephase.cpp \
    particlephase/dynamicparticlephase.cpp \
    particlephase/particlephase.cpp \
    particlephase/pbmbccaparticlephase.cpp \
    particlephase/pbmdlcaparticlephase.cpp \
    particlephase/randomwalkparticlephase.cpp \
    particlephase/rigidbodylangevinparticlephase.cpp \
    dynamicparticle.cpp \
    alglib/src/alglibinternal.cpp \
    alglib/src/alglibmisc.cpp \
    alglib/src/ap.cpp \
    alglib/src/dataanalysis.cpp \
    alglib/src/diffequations.cpp \
    alglib/src/fasttransforms.cpp \
    alglib/src/integration.cpp \
    alglib/src/interpolation.cpp \
    alglib/src/linalg.cpp \
    alglib/src/optimization.cpp \
    alglib/src/solvers.cpp \
    alglib/src/specialfunctions.cpp \
    alglib/src/statistics.cpp \
    spline.cpp \
    splinelinear.cpp \
    particlephase/pbmfractalparticlephase.cpp \
    particlephase/pbmparticlephase.cpp \
    grid.cpp \
    cg_simulation_dataXML.cpp \
    XML_config_file.cpp \
    linked_langevin_simulation.cpp \
    linked_moments_simulation.cpp \
    linked_pbm_simulation.cpp \
    linked_simulation.cpp \
    meso_simulation.cpp \
    stand_alone_langevin_simulation.cpp \
    stand_alone_moments_simulation.cpp \
    stand_alone_pbm_simulation.cpp \
    stand_alone_simulation.cpp \
    sectional/sectional.cpp \
    dynamic_timestep.cpp \
    xml_templates.cpp \
    linked_moments_temp_simulation.cpp \
    linked_pbm_temp_simulation.cpp

HEADERS += \
    nanodome.h \
    species.h \
    particle.h \
    mesoobject.h \
    edge.h \
    utilities.h \
    virtualbond.h \
    pointbond.h \
    fixedbond.h \
    dynamicparticle.h \
    pointparticle.h \
    cnt.h \
    gasphase.h \
    nucleationtheory.h \
    particlephase.h \
    point.h \
    dynamicpoint.h \
    ndm_random.h \
    objectcounter.h \
    bond.h \
    clock.h \
    particlebond.h \
    collisionkernel.h \
    collisionalobject.h \
    physicalobject.h \
    coagulationkernel.h \
    binarykernel.h \
    fmbinarykernel.h \
    majorantkernel.h \
    fmmajorantkernel.h \
    tetra_bond.h \
    tetra_constants.h \
    tetra_constrainer.h \
    tetra_constraint.h \
    tetra_face.h \
    tetra_link.h \
    tetra_tetrahedron.h \
    tetra_vertex.h \
    tinyxml2/tinyxml2.h \
    i_o.h \
    JSONFile.h \
    streamline.h \
    streamlinefileJSON.h \
    streamlinefileXML.h \
    JSON/JSON.h \
    JSON/JSONValue.h \
    XMLfile.h \
    aggregate/aggregate.h \
    aggregate/ellipsoidaggregate.h \
    aggregate/fractalaggregate.h \
    aggregate/fullrigidbodyaggregate.h \
    aggregate/pbmaggregate.h \
    aggregate/rattleaggregate.h \
    aggregate/rigidbodyaggregate.h \
    aggregate/spatialaggregate.h \
    gasphase/gasphase.h \
    gasphase/gasphasecv.h \
    gasphase/gasphaselink.h \
    moments/momentmodel.h \
    moments/momentmodelfriedlander.h \
    moments/momentmodelpratsinis.h \
    particlephase/cgmdparticlephase.h \
    particlephase/constrainedlangevinparticlephase.h \
    particlephase/dynamicparticlephase.h \
    particlephase/particlephase.h \
    particlephase/pbmbccaparticlephase.h \
    particlephase/pbmdlcaparticlephase.h \
    particlephase/randomwalkparticlephase.h \
    particlephase/rigidbodylangevinparticlephase.h \
    alglib/src/alglibinternal.h \
    alglib/src/alglibmisc.h \
    alglib/src/ap.h \
    alglib/src/dataanalysis.h \
    alglib/src/diffequations.h \
    alglib/src/fasttransforms.h \
    alglib/src/integration.h \
    alglib/src/interpolation.h \
    alglib/src/linalg.h \
    alglib/src/optimization.h \
    alglib/src/solvers.h \
    alglib/src/specialfunctions.h \
    alglib/src/statistics.h \
    alglib/src/stdafx.h \
    spline.h \
    splinelinear.h \
    particlephase/pbmfractalparticlephase.h \
    particlephase/pbmparticlephase.h \
    grid.h \
    cg_simulation_dataXML.h \
    XML_config_file.h \
    linked_langevin_simulation.h \
    linked_moments_simulation.h \
    linked_pbm_simulation.h \
    linked_simulation.h \
    stand_alone_langevin_simulation.h \
    stand_alone_moments_simulation.h \
    stand_alone_pbm_simulation.h \
    stand_alone_simulation.h \
    meso_simulation.h \
    sectional/sectional.h \
    dynamic_timestep.h \
    xml_templates.h \
    linked_moments_temp_simulation.h \
    linked_pbm_temp_simulation.h

