from importlib import import_module, reload

import pytest

from ga4gh.vrs import config


@pytest.fixture
def configured_modules(monkeypatch, request):
    """Reload import-time defaults, restoring them after each environment case."""
    modules = [
        config,
        import_module("ga4gh.vrs.normalize"),
        import_module("ga4gh.vrs.extras.translator"),
        import_module("ga4gh.vrs.extras.annotator.vcf"),
    ]
    original_namespaces = [vars(module).copy() for module in modules]
    try:
        with monkeypatch.context() as env:
            env.delenv("GA4GH_VRS_RLE_SEQ_LIMIT", raising=False)
            if request.param is not None:
                env.setenv("GA4GH_VRS_RLE_SEQ_LIMIT", request.param)
            for module in modules:
                reload(module)
            yield modules[1:]
    finally:
        # Preserve class identities already imported by other test modules.
        for module, namespace in zip(modules, original_namespaces, strict=True):
            vars(module).clear()
            vars(module).update(namespace)
