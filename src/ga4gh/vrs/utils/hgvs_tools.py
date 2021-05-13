import hgvs


class HgvsTools():

    """ A utility class that provides access to hgvs components that have a UTA database dependency"""
    def __init__(self, uta_conn):

        self.uta_conn = uta_conn
        self.normalizer = hgvs.normalizer.Normalizer(self.uta_conn)
        self.parser = hgvs.parser.Parser()

    def close(self):
        # TODO These should only be closed if they are owned by this instance
        uta_conn = self.uta_conn
        self.normalizer = None
        if uta_conn is not None:
            uta_conn.close()

    def parse(self, hgvs_str):
        return self.parser.parse_hgvs_variant(hgvs_str)

    def normalize(self, hgvs):
        norm_hgvs = self.normalizer.normalize(hgvs)

        return norm_hgvs