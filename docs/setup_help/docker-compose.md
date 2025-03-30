The compose.yaml file in this repo can be used by Docker Compose. It has also been tested to work with Podman Compose.

See full specification:
https://github.com/compose-spec/compose-spec/blob/main/00-overview.md


The default docker compose in this repo assumes you have not downloaded the UTA or SeqRepo databases out of band, but either ran the compose from scratch which runs the `biocommons/seqrepo` or `biocommons/uta` containers in a way that populates their databases if they don't exist, or you ran those containers on your own previously, such that the volumes `seqrepo_vol` and `uta_vol` are already populated.

If you already have a SeqRepo database directory on your local filesystem, you can mount it directly as a volume, or can copy it into a docker volume (recommended if your local disk space is not a concern).


For an example  of mounting a local seqrepo dir directly to the seqrepo-rest-service container see [seqrepo-mount-local.compose.yaml](./docker-compose-examples/seqrepo-mount-local.compose.yaml)


For an example  of populating a docker volume with a local seqrepo dir see [seqrepo-copy-local.compose.yaml](./docker-compose-examples/seqrepo-copy-local.compose.yaml)


In both of the above examples, a volume called `uta_dl_cache` is used. UTA downloads a postgres dump archive if its database is not populated yet, and stores the archive in `/tmp`. If we make `/tmp` a persistent volume, we can destroy the UTA container and the `uta_vol` and recreate it from scratch using the compressed archive in `uta_dl_cache`. The latest UTA compressed postgres dump for `biocommons/uta:uta_20241220` is 344MB, while the uncompressed postgres database created from it is 6GB.

To run one of the compose files that uses a local seqrepo, with a seqrepo root dir in `~/.local/share/seqrepo/`, run

```
$ SEQREPO_ROOT_DIR=$HOME/.local/share/seqrepo/ docker-compose -f $(pwd)/docs/setup_help/docker-compose-examples/seqrepo-copy-local.compose.yaml up
```
