from .bundlemanager import BundleManager



class DemoApp(BundleManager):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.annotations = {}

    
