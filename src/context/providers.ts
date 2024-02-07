import { ContextProvider } from "@exabyte-io/code.js/dist/context";
import { ContextProviderConfig } from "@exabyte-io/code.js/dist/context/provider";

export class ExecutableContextProvider extends ContextProvider {
    constructor(config: ContextProviderConfig) {
        super({
            ...config,
            domain: "executable",
        });
    }
}
