#Notes from Wes Goar on setup of vrs-python on a m1 mac
- date: 20220328
- macOS version: 12.3
- python version: 3.9.11

##Notes
- I installed seqrepo before trying this setup.
- I had already installed xcode and homebrew on my machine.
- This installation assumes you are using zsh. Please change the commands to fit your environment of choice.

###Steps
1. brew install openssl
2. add the following statement in your .zshrc:
```shell
export PATH="/opt/homebrew/opt/openssl@1.1/bin:$PATH"
```

3. Run brew info openssl to ascertain the correct environment path to put into the export statements (add the following statements, but follow your info, into your .zshenv:)
```shell
export LDFLAGS="-L/opt/homebrew/opt/openssl@1.1/lib"
export CPPFLAGS="-I/opt/homebrew/opt/openssl@3/include"
export PKG_CONFIG_PATH="/opt/homebrew/opt/openssl@1.1/lib/pkgconfig"
```

4. source ~/.zshrc
5. source ~/.zshenv
6. brew install libpq
7. add the following statement in your .zshrc:
```shell
export PATH="/opt/homebrew/opt/libpq/bin:$PATH"
```
8. Run brew info libpq to ascertain the correct environment path to append to the appropriate flags (add the following statements, but follow your info, into your .zshenv:)

```shell
export LDFLAGS="-L/opt/homebrew/opt/openssl@1.1/lib -L/opt/homebrew/opt/libpq/lib"
export CPPFLAGS="-I/opt/homebrew/opt/openssl@3/include -I/opt/homebrew/opt/libpq/include"
export PKG_CONFIG_PATH="/opt/homebrew/opt/openssl@1.1/lib/pkgconfig:/opt/homebrew/opt/libpq/lib/pkgconfig"
```

9. source ~/.zshenv
10. brew install postgres
11. add the following statement in your .zshrc: export PATH="/opt/homebrew/opt/postgresql@14/bin:$PATH"
    ###Make sure that you update the @14 with your own version
12. source ~/.zshrc
13. Install UTA (/docs/setup_help/uta_installation.md)
14. Run the make devready command:
    1. `make devready`
15. Run the make test command:
    1. `make test`