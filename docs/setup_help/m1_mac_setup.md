#Notes from Wes Goar on setup of vrs-python on a m1 mac
- @author Wes Goar
- date: 20220328
- macOS version: 12.3
- python version: 3.9.11

##Notes
- I installed seqrepo before trying this setup.
- I had already installed xcode and homebrew on my machine.

##Steps
1. `brew install openssl`
2. `brew install libpq`
3. `brew install postgres`
4. Create .zshrc and .zshenv files if they do not exist already
5. Add the following to your .zshrc file (**_NOTE: Your path/version numbers may differ from what is shown here. Use the ones on your machine._** You may get this info with: `brew info thinghere`): 

    ```
    export PATH="/opt/homebrew/opt/openssl@3/bin:$PATH"
    export PATH="/opt/homebrew/opt/libpq/bin:$PATH"
    export PATH="/opt/homebrew/opt/postgresql@14/bin:$PATH"
   ```
6. Add the following to your .zshenv file (make sure to update the correct version/path numbers here as well)
    ```
    export LDFLAGS="-L/opt/homebrew/opt/openssl@3/lib -L/opt/homebrew/opt/libpq/lib"
    export CPPFLAGS="-I/opt/homebrew/opt/openssl@3/include -I/opt/homebrew/opt/libpq/include"
    export PKG_CONFIG_PATH="/opt/homebrew/opt/openssl@3/lib/pkgconfig:/opt/homebrew/opt/libpq/lib/pkgconfig"
      ```

7. [Install UTA](../uta_installation.md)
8. Return to the vrs-python installation steps. You're almost done! 


