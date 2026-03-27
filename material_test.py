# To run:
# pytest material_test.py

import pytest
from material import concrete_stress_strain, steel_rebar_stress_strain

def test_concrete_stress_strain():
    assert concrete_stress_strain(fc__=28, gamma=0.15, eps_0=0.002, eps_u=0.003, eps_c=0.000) == pytest.approx(0, abs=0.01)
    assert concrete_stress_strain(fc__=28, gamma=0.15, eps_0=0.002, eps_u=0.003, eps_c=0.0005) == pytest.approx(12.25, abs=0.01)
    assert concrete_stress_strain(fc__=28, gamma=0.15, eps_0=0.002, eps_u=0.003, eps_c=0.001) == pytest.approx(21, abs=0.01)
    assert concrete_stress_strain(fc__=28, gamma=0.15, eps_0=0.002, eps_u=0.003, eps_c=0.0015) == pytest.approx(26.25, abs=0.01)
    assert concrete_stress_strain(fc__=28, gamma=0.15, eps_0=0.002, eps_u=0.003, eps_c=0.002) == pytest.approx(28, abs=0.01)
    assert concrete_stress_strain(fc__=28, gamma=0.15, eps_0=0.002, eps_u=0.003, eps_c=0.0025) == pytest.approx(25.9, abs=0.01)
    assert concrete_stress_strain(fc__=28, gamma=0.15, eps_0=0.002, eps_u=0.003, eps_c=0.003) == pytest.approx(23.8, abs=0.01)


def test_steel_rebar_stress_strain():
    assert steel_rebar_stress_strain(Es=200000, fy=280, fsu=400, eps_sh=0.015, eps_su=0.05, eps_s=0.000) == pytest.approx(0, abs=0.01)
    assert steel_rebar_stress_strain(Es=200000, fy=280, fsu=400, eps_sh=0.015, eps_su=0.05, eps_s=0.0005) == pytest.approx(100, abs=0.01)
    assert steel_rebar_stress_strain(Es=200000, fy=280, fsu=400, eps_sh=0.015, eps_su=0.05, eps_s=0.001) == pytest.approx(200, abs=0.01)
    assert steel_rebar_stress_strain(Es=200000, fy=280, fsu=400, eps_sh=0.015, eps_su=0.05, eps_s=0.0015) == pytest.approx(280, abs=0.01)
    assert steel_rebar_stress_strain(Es=200000, fy=280, fsu=400, eps_sh=0.015, eps_su=0.05, eps_s=0.002) == pytest.approx(280, abs=0.01)
    assert steel_rebar_stress_strain(Es=200000, fy=280, fsu=400, eps_sh=0.015, eps_su=0.05, eps_s=0.015) == pytest.approx(280, abs=0.01)
    assert steel_rebar_stress_strain(Es=200000, fy=280, fsu=400, eps_sh=0.015, eps_su=0.05, eps_s=0.02) == pytest.approx(323.33, abs=0.01)
    assert steel_rebar_stress_strain(Es=200000, fy=280, fsu=400, eps_sh=0.015, eps_su=0.05, eps_s=0.025) == pytest.approx(352.90, abs=0.01)
    assert steel_rebar_stress_strain(Es=200000, fy=280, fsu=400, eps_sh=0.015, eps_su=0.05, eps_s=0.03) == pytest.approx(372.97, abs=0.01)
    assert steel_rebar_stress_strain(Es=200000, fy=280, fsu=400, eps_sh=0.015, eps_su=0.05, eps_s=0.04) == pytest.approx(394.4, abs=0.01)
    assert steel_rebar_stress_strain(Es=200000, fy=280, fsu=400, eps_sh=0.015, eps_su=0.05, eps_s=0.05) == pytest.approx(400, abs=0.01)
