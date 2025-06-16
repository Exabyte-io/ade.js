import { ContextProvider } from "@mat3ra/code/dist/js/context";
import type { ContextProviderConfig } from "@mat3ra/code/dist/js/context/provider";

export class ExecutableContextProvider extends ContextProvider {
    constructor(config: ContextProviderConfig) {
        super({
            ...config,
            domain: "executable",
        });
    }
}
